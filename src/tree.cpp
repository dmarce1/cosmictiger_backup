#include <cosmictiger/defs.hpp>
#include <cosmictiger/util.hpp>
#include <cosmictiger/tree.hpp>
#include <cosmictiger/tree_dir.hpp>

#include <atomic>
#include <stack>
#include <thread>

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<tree>, tree);

using build_tree_dir_action_type = tree::build_tree_dir_action;
using destroy_action_type = tree::destroy_action;
using drift_action_type = tree::drift_action;
using find_home_action_type = tree::find_home_action;
using grow_action_type = tree::grow_action;
using load_balance_action_type = tree::load_balance_action;
using prune_action_type = tree::prune_action;
using get_parts_action_type = tree::get_parts_action;
using get_ptr_action_type = tree::get_ptr_action;
using migrate_action_type = tree::migrate_action;
using verify_action_type = tree::verify_action;

HPX_REGISTER_ACTION (build_tree_dir_action_type);
HPX_REGISTER_ACTION (destroy_action_type);
HPX_REGISTER_ACTION (drift_action_type);
HPX_REGISTER_ACTION (find_home_action_type);
HPX_REGISTER_ACTION (grow_action_type);
HPX_REGISTER_ACTION (load_balance_action_type);
HPX_REGISTER_ACTION (prune_action_type);
HPX_REGISTER_ACTION (get_parts_action_type);
HPX_REGISTER_ACTION (get_ptr_action_type);
HPX_REGISTER_ACTION (migrate_action_type);
HPX_REGISTER_ACTION (verify_action_type);

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

tree::tree(box_id_type id) {
	tptr = new tree_mems;
	tptr->boxid = id;
	tptr->leaf = true;
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
	const int my_level = msb(tptr->boxid) - 1;
	if (my_level == min_level) {
		tree_dir dir;
		dir.add_tree_client(self, tptr->boxid);
		return dir;
	} else {
		auto futl = hpx::async < build_tree_dir_action > (tptr->children[0].get_id(), tptr->children[0]);
		auto futr = hpx::async < build_tree_dir_action > (tptr->children[1].get_id(), tptr->children[1]);
		tree_dir left = futl.get();
		return left.merge(futr.get());
	}
}

void tree::create_children() {
	auto futl = hpx::new_ < tree > (hpx::find_here(), (tptr->boxid << box_id_type(1)));
	auto futr = hpx::new_ < tree > (hpx::find_here(), (tptr->boxid << box_id_type(1)) + box_id_type(1));
	tptr->leaf = false;
	auto id_left = futl.get();
	auto id_right = futr.get();
	auto fptrl = hpx::async < get_ptr_action > (id_left);
	auto fptrr = hpx::async < get_ptr_action > (id_right);
	tptr->children[0] = tree_client(std::move(id_left), fptrl.get());
	tptr->children[1] = tree_client(std::move(id_right), fptrr.get());
}

int tree::drift(int stack_cnt, int step, tree_client parent, tree_client self, float dt) {
	assert(tptr->boxid);
	tptr->parent = parent;
	if (tptr->leaf) {
		std::unique_lock<mutex_type> lock(tptr->mtx);
		bucket exit_parts;
		const auto box = box_id_to_range(tptr->boxid);
		auto i = tptr->parts.begin();
		while (i != tptr->parts.end()) {
			if (step % 2 == i->step) {
				auto x = pos_to_double(i->x);
				const auto v = i->v;
				x += v * dt;
				i->x = double_to_pos(x);
				i->step++;
				if (!in_range(pos_to_double(i->x), box)) {
					exit_parts.insert(*i);
					i = tptr->parts.remove(i);
				} else {
					i++;
				}
			} else {
				i++;
			}
		}
		lock.unlock();
		if (exit_parts.size()) {
			parent.find_home(stack_cnt, std::move(exit_parts)).get();
		}
	} else {
		auto futl = tptr->children[0].drift(stack_cnt, step, self, tptr->children[0], dt);
		auto futr = tptr->children[1].drift(stack_cnt, step, self, tptr->children[1], dt);
		futl.get();
		futr.get();
	}
	return 0;
}

int tree::find_home(int stack_cnt, bucket parts) {

	const auto box = box_id_to_range(tptr->boxid);
	std::unique_lock<mutex_type> lock(tptr->mtx);
	if (tptr->leaf) {
		while (parts.size()) {
			assert(in_range(pos_to_double(parts.front().x), box));
			tptr->parts.insert(parts.front());
			parts.remove(parts.begin());
		}
		lock.unlock();
	} else {
		bucket p_parts;
		bucket l_parts;
		bucket r_parts;
		const auto boxl = box_id_to_range(tptr->boxid << box_id_type(1));
		const auto boxr = box_id_to_range((tptr->boxid << box_id_type(1)) + box_id_type(1));
		while (parts.size()) {
			auto &p = parts.front();
			const auto x = pos_to_double(p.x);
			if (in_range(x, boxl)) {
				l_parts.insert(p);
			} else if (in_range(x, boxr)) {
				r_parts.insert(p);
			} else {
				p_parts.insert(p);
			}
			parts.remove(parts.begin());
		}
		lock.unlock();
		hpx::future<int> pfut;
		hpx::future<int> lfut;
		hpx::future<int> rfut;
		if (p_parts.size()) {
			assert(tptr->parent != tree_client());
			pfut = tptr->parent.find_home(stack_cnt, std::move(p_parts));
		} else {
			pfut = hpx::make_ready_future(0);
		}
		if (l_parts.size()) {
			assert(tptr->children[0] != tree_client());
			lfut = tptr->children[0].find_home(stack_cnt, std::move(l_parts));
		} else {
			lfut = hpx::make_ready_future(0);
		}
		if (r_parts.size()) {
			assert(tptr->children[1] != tree_client());
			rfut = tptr->children[1].find_home(stack_cnt, std::move(r_parts));
		} else {
			rfut = hpx::make_ready_future(0);
		}
		pfut.get();
		lfut.get();
		rfut.get();
	}
	return 0;
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
	const int my_level = msb(tptr->boxid) - 1;
	if (tptr->leaf) {
		if (my_level < min_level || parts.size() + tptr->parts.size() > opts.bucket_size) {
			create_children();
		} else {
			while (parts.size()) {
				tptr->parts.insert(parts.front());
				parts.remove(parts.begin());
			}
		}
	}
	if (!tptr->leaf) {
		bucket parts_left;
		bucket parts_right;
		while (tptr->parts.size()) {
			parts.insert(tptr->parts.front());
			tptr->parts.remove(tptr->parts.begin());
		}
		range box = box_id_to_range(tptr->boxid);
		const int dim = range_max_dim(box);
		const double xmid = 0.5 * (box.max[dim] + box.min[dim]);
		while (parts.size()) {
			const auto &p = parts.front();
			if (double(p.x[dim]) < xmid) {
				parts_left.insert(p);
			} else {
				parts_right.insert(p);
			}
			parts.remove(parts.begin());
		}
//		if (!tptr->children[0].local()) {
//			printf("Sending %i left\n", parts_left.size());
//		}
//		if (!tptr->children[1].local()) {
//			printf("Sending %i right\n", parts_right.size());
//		}
		auto futl = tptr->children[0].grow(stack_cnt, std::move(parts_left));
		auto futr = tptr->children[1].grow(stack_cnt, std::move(parts_right));
		tptr->child_cnt[0] = futl.get();
		tptr->child_cnt[1] = futr.get();
	}
//	printf("%x\n", tptr->boxid);
	return size();
}

bucket tree::get_parts() {
	bucket parts;
	while (tptr->parts.size()) {
		parts.insert(tptr->parts.front());
		tptr->parts.remove(tptr->parts.begin());
	}
	if (!tptr->leaf) {
		auto futl = tptr->children[0].get_parts();
		auto futr = tptr->children[1].get_parts();
		auto tmp = futl.get();
		while (tmp.size()) {
			parts.insert(tmp.front());
			tmp.remove(tmp.begin());
		}
		tmp = futr.get();
		while (tmp.size()) {
			parts.insert(tmp.front());
			tmp.remove(tmp.begin());
		}
	}
	return parts;
}

int tree::load_balance(int stack_cnt, std::uint64_t index) {
	if (!tptr->leaf) {
		const auto &localities = hpx_localities();
		const int min_level = msb(localities.size()) - 1;
		std::uint64_t il, ir;
		const int child_level = msb(tptr->boxid);
		if (child_level > min_level) {
			il = index * std::uint64_t(localities.size()) / opts.problem_size;
			ir = (index + tptr->child_cnt[0]) * std::uint64_t(localities.size()) / opts.problem_size;
			il = hpx::get_locality_id();
			ir = hpx::get_locality_id();
		} else {
			const auto total_nodes = (1 << child_level);
			const auto boxl = tptr->boxid << box_id_type(1);
			const auto boxr = (tptr->boxid << box_id_type(1)) + box_id_type(1);
			il = (boxl - total_nodes) * localities.size() / total_nodes;
			ir = (boxr - total_nodes) * localities.size() / total_nodes;
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
			while (tmpl.size()) {
				tptr->parts.insert(tmpl.front());
				tmpl.remove(tmpl.begin());
			}
			while (tmpr.size()) {
				tptr->parts.insert(tmpr.front());
				tmpr.remove(tmpr.begin());
			}
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
	const int my_level = msb(tptr->boxid) - 1;
	int rc = 0;
	if (tptr->leaf) {
		if (size() > opts.bucket_size) {
			rc |= TREE_OVERFLOW;
		}
		const auto box = box_id_to_range(tptr->boxid);
		for (auto iter = tptr->parts.begin(); iter != tptr->parts.end(); iter++) {
			const auto this_x = pos_to_double(iter->x);
			if (!in_range(this_x, box)) {
				rc |= TREE_INVALID;
			}
		}
	} else {
		if (my_level > min_level && size() <= opts.bucket_size) {
			rc |= TREE_UNDERFLOW;
		}
		auto futl = tptr->children[0].verify(stack_cnt);
		auto futr = tptr->children[1].verify(stack_cnt);
		rc |= futl.get();
		rc |= futr.get();
	}
	return rc;
}

