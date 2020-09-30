#include <cosmictiger/defs.hpp>
#include <cosmictiger/util.hpp>
#include <cosmictiger/tree.hpp>
#include <cosmictiger/tree_dir.hpp>

#include <atomic>
#include <stack>
#include <thread>

HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);

std::string tree_verification_error(int rc) {
	std::string error = "";
	if (rc & TREE_OVERFLOW) {
		error += "TREE_OVERFLOW ";
	}
	if (rc & TREE_UNDERFLOW) {
		error += "TREE_UNDERFLOW ";
	}
	if (rc & TREE_INVALID) {
		error += "TREE_INVALID ";
	}
	return error;
}

tree::tree() {
	tptr = new tree_mems;
}

tree::tree(const range &box, int level) {
	tptr = new tree_mems;
	tptr->box = box;
	tptr->leaf = true;
	tptr->level = level;
}

tree::tree(const tree &other) {
	tptr = new tree_mems;
	*tptr = *(other.tptr);
}

tree::~tree() {
	delete tptr;
}

tree_dir tree::build_tree_dir(tree_client self) const {
	const int min_level = msb(hpx_localities().size()) - 1;
	if (tptr->level == min_level) {
		tree_dir dir;
		dir.add_tree_client(self, tptr->box);
		return dir;
	} else {
		auto futl = hpx::async < build_tree_dir_action > (tptr->children[0].get_id(), tptr->children[0]);
		auto futr = hpx::async < build_tree_dir_action > (tptr->children[1].get_id(), tptr->children[1]);
		tree_dir left = futl.get();
		return left.merge(futr.get());
	}
}

void tree::create_children() {
	auto boxl = tptr->box;
	auto boxr = tptr->box;
	const auto dim = tptr->level % NDIM;
	boxl.max[dim] = boxr.min[dim] = 0.5 * (tptr->box.min[dim] + tptr->box.max[dim]);
	auto futl = hpx::new_ < tree > (hpx::find_here(), boxl, tptr->level + 1);
	auto futr = hpx::new_ < tree > (hpx::find_here(), boxr, tptr->level + 1);
	tptr->leaf = false;
	auto id_left = futl.get();
	auto id_right = futr.get();
	auto fptrl = hpx::async < get_ptr_action > (id_left);
	auto fptrr = hpx::async < get_ptr_action > (id_right);
	tptr->children[0] = tree_client(std::move(id_left), fptrl.get());
	tptr->children[1] = tree_client(std::move(id_right), fptrr.get());
}

int tree::destroy(int stack_cnt) {
	if (!tptr->leaf) {
		auto futl = tptr->children[0].destroy(stack_cnt);
		auto futr = tptr->children[1].destroy(stack_cnt);
		futl.get();
		tptr->children[0] = tree_client();
		futr.get();
		tptr->children[1] = tree_client();
	}
	return 0;
}

std::uint64_t tree::get_ptr() {
	return reinterpret_cast<std::uint64_t>(this);
}

std::uint64_t tree::grow(int stack_cnt, bucket &&parts) {
	const int min_level = msb(hpx_localities().size()) - 1;
	if (tptr->leaf) {
		if (tptr->level < min_level || parts.size() + tptr->parts.size() > opts.bucket_size) {
			create_children();
		} else {
			if (tptr->parts.size() == 0) {
				tptr->parts = std::move(parts);
			} else {
				for (auto i = parts.begin(); i != parts.end(); i++) {
					tptr->parts.insert(*i);
				}
			}
		}
	}
	if (!tptr->leaf) {
		while (tptr->parts.size()) {
			parts.insert(tptr->parts.front());
			tptr->parts.remove(tptr->parts.begin());
		}
		bucket parts_left = std::move(parts);
		bucket parts_right;
		const int dim = range_max_dim(tptr->box);
		const double xmid = 0.5 * (tptr->box.max[dim] + tptr->box.min[dim]);
		auto i = parts_left.begin();
		while (i != parts_left.end()) {
			if (double(i->x[dim]) >= xmid) {
				parts_right.insert(*i);
				i = parts_left.remove(i);
			} else {
				i++;
			}
		}
		auto futl = tptr->children[0].grow(stack_cnt, std::move(parts_left));
		auto futr = tptr->children[1].grow(stack_cnt, std::move(parts_right));
		tptr->child_cnt[0] = futl.get();
		tptr->child_cnt[1] = futr.get();
	}
	return size();
}

bucket tree::get_parts() {
	return std::move(tptr->parts);
}

int tree::load_balance(int stack_cnt, std::uint64_t index) {
	if (!tptr->leaf) {
		const auto &localities = hpx_localities();
		const int min_level = msb(localities.size()) - 1;
		std::uint64_t il, ir;
		const int child_level = tptr->level + 1;
		if (child_level > min_level) {
			il = index * std::uint64_t(localities.size()) / opts.problem_size;
			ir = (index + tptr->child_cnt[0]) * std::uint64_t(localities.size()) / opts.problem_size;
		} else {
			const auto total_nodes = (1 << child_level);
			const auto dim = tptr->level % NDIM;
			const auto xl = 0.75 * tptr->box.min[dim] + 0.25 * tptr->box.max[dim];
			const auto xr = 0.25 * tptr->box.min[dim] + 0.75 * tptr->box.max[dim];
			const auto x0 = (1 << (tptr->level / NDIM));
			il = std::uint64_t(xl * x0) * localities.size() / total_nodes;
			ir = std::uint64_t(xr * x0) * localities.size() / total_nodes;
		}
		assert(il >= 0);
		assert(il < localities.size());
		assert(ir >= 0);
		assert(ir < localities.size());
		hpx::future<tree_client> newl;
		hpx::future<tree_client> newr;
		if (il != hpx::get_locality_id()) {
			newl = tptr->children[0].migrate(localities[il]);
		}
		if (ir != hpx::get_locality_id()) {
			newr = tptr->children[1].migrate(localities[ir]);
		}
		if (newl.valid()) {
			tptr->children[0] = newl.get();
		}
		if (newr.valid()) {
			tptr->children[1] = newr.get();
		}
		auto futl = tptr->children[0].load_balance(stack_cnt, index);
		auto futr = tptr->children[1].load_balance(stack_cnt, index + tptr->child_cnt[0]);
		futl.get();
		futr.get();
	}
	return 0;
}

tree_client tree::migrate(hpx::id_type locality) {
	auto new_ptr = hpx::new_ < tree > (locality, *this);
	auto new_id = new_ptr.get();
	return tree_client(new_id, hpx::async < get_ptr_action > (new_id).get());

}

std::uint64_t tree::prune(int stack_cnt) {
	tptr->parent = tree_client();
	if (!tptr->leaf) {
		auto futl = tptr->children[0].prune(stack_cnt);
		auto futr = tptr->children[1].prune(stack_cnt);
		tptr->child_cnt[0] = futl.get();
		tptr->child_cnt[1] = futr.get();
		if (size() <= opts.bucket_size) {
			auto futl = tptr->children[0].get_parts();
			auto futr = tptr->children[1].get_parts();
			auto tmpl = futl.get();
			tptr->children[0] = tree_client();
			auto tmpr = futr.get();
			tptr->children[1] = tree_client();
			auto sz = tmpr.size() + tmpl.size();
			for (auto i = tmpl.begin(); i != tmpl.end(); i++) {
				tmpr.insert(*i);
			}
			tptr->parts = std::move(tmpr);
			tptr->leaf = true;
		}
	} else if (tptr->parts.size() > opts.bucket_size) {
		grow(stack_cnt + 1, bucket());
	}
	return size();
}

std::size_t tree::size() const {
	if (tptr->leaf) {
		return tptr->parts.size();
	} else {
		return tptr->child_cnt[0] + tptr->child_cnt[1];
	}
}

int tree::verify(int stack_cnt) const {
	const int min_level = msb(hpx_localities().size()) - 1;
	int rc = 0;
	if (tptr->leaf) {
		if (size() > opts.bucket_size) {
			rc |= TREE_OVERFLOW;
		}
		for (auto iter = tptr->parts.begin(); iter != tptr->parts.end(); iter++) {
			const auto this_x = pos_to_double(iter->x);
			if (!in_range(this_x, tptr->box)) {
				rc |= TREE_INVALID;
			}
		}
	} else {
		if (tptr->level > min_level && size() <= opts.bucket_size) {
			rc |= TREE_UNDERFLOW;
		}
		auto futl = tptr->children[0].verify(stack_cnt);
		auto futr = tptr->children[1].verify(stack_cnt);
		rc |= futl.get();
		rc |= futr.get();
	}
	return rc;
}

