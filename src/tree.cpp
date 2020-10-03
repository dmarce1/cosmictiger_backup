#include <cosmictiger/defs.hpp>
#include <cosmictiger/util.hpp>
#include <cosmictiger/tree.hpp>
#include <cosmictiger/tree_dir.hpp>

#include <atomic>
#include <stack>
#include <thread>

#ifdef HPX_LITE
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<tree>, tree);
using build_tree_dir_type = tree::build_tree_dir_action;
using destroy_type = tree::destroy_action;
using drift_type = tree::drift_action;
using find_home_parent_type = tree::find_home_parent_action;
using find_home_child_type = tree::find_home_child_action;
using grow_type = tree::grow_action;
using load_balance_type = tree::load_balance_action;
using prune_type = tree::prune_action;
using verify_type = tree::verify_action;
using get_child_checks_type = tree::get_child_checks_action;
using get_ptr_type = tree::get_ptr_action;
using get_parts_type = tree::get_parts_action;
using migrate_type = tree::migrate_action;
HPX_REGISTER_ACTION(build_tree_dir_type);
HPX_REGISTER_ACTION(destroy_type);
HPX_REGISTER_ACTION(drift_type);
HPX_REGISTER_ACTION(find_home_parent_type);
HPX_REGISTER_ACTION(find_home_child_type);
HPX_REGISTER_ACTION(grow_type);
HPX_REGISTER_ACTION(load_balance_type);
HPX_REGISTER_ACTION(prune_type);
HPX_REGISTER_ACTION(verify_type);
HPX_REGISTER_ACTION(get_child_checks_type);
HPX_REGISTER_ACTION(get_ptr_type);
HPX_REGISTER_ACTION(get_parts_type);
HPX_REGISTER_ACTION(migrate_type);

#else
HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);
#endif

static std::shared_ptr<tree_dir> directory;

void tree_broadcast_directory(const tree_dir &dir);
void tree_cleanup();
void tree_complete_drift();

HPX_PLAIN_ACTION (tree_broadcast_directory);
HPX_PLAIN_ACTION (tree_complete_drift);
HPX_PLAIN_ACTION (tree_cleanup);

void tree_insert_parts(bucket&& parts) {
	directory->find_home(0,std::move(parts));
	directory->retire_futures();
}

void tree_broadcast_directory(const tree_dir &dir) {
	const auto &localities = hpx_localities();
	const int il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const int ir = ((hpx::get_locality_id() + 1) << 1);
	std::vector<hpx::future<void>> futs;
	if (il < localities.size()) {
		futs.push_back(hpx::async < tree_broadcast_directory_action > (localities[il], dir));
	}
	if (ir < localities.size()) {
		futs.push_back(hpx::async < tree_broadcast_directory_action > (localities[ir], dir));
	}
//	printf( "Setting dir %i\n", dir.size());
	directory = std::make_shared < tree_dir > (dir);
	hpx::wait_all(futs.begin(), futs.end());
}

void tree_complete_drift() {
	const auto &localities = hpx_localities();
	const int il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const int ir = ((hpx::get_locality_id() + 1) << 1);
	std::vector<hpx::future<void>> futs;
	if (il < localities.size()) {
		futs.push_back(hpx::async < tree_complete_drift_action > (localities[il]));
	}
	if (ir < localities.size()) {
		futs.push_back(hpx::async < tree_complete_drift_action > (localities[ir]));
	}
	directory->retire_futures();
	hpx::wait_all(futs.begin(), futs.end());
}

void tree_cleanup() {
	const auto &localities = hpx_localities();
	const int il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const int ir = ((hpx::get_locality_id() + 1) << 1);
	std::vector<hpx::future<void>> futs;
	if (il < localities.size()) {
		futs.push_back(hpx::async < tree_cleanup_action > (localities[il]));
	}
	if (ir < localities.size()) {
		futs.push_back(hpx::async < tree_cleanup_action > (localities[ir]));
	}
	directory = nullptr;
	hpx::wait_all(futs.begin(), futs.end());
}

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

tree::tree(box_id_type id, const range &box, int level) {
	tptr = new tree_mems;
	tptr->box = box;
	tptr->leaf = true;
	tptr->level = level;
	tptr->id = id;
}

tree::tree(const tree &other) {
	tptr = new tree_mems;
	*tptr = *(other.tptr);
}

tree::~tree() {
	delete tptr;
}

tree_dir tree::build_tree_dir(tree_client self) const {
	const int min_level = bits_to_level(hpx_localities().size());
	tree_dir dir;
	if (tptr->level == min_level) {
		dir.add_tree_client(self, tptr->box);
	} else {
		auto futl = hpx::async < build_tree_dir_action > (tptr->left_child.get_id(), tptr->left_child);
		auto futr = hpx::async < build_tree_dir_action > (tptr->right_child.get_id(), tptr->right_child);
		dir = futl.get();
//		printf( "Merging\n");
		dir.merge(futr.get());
	}
	if (tptr->level == 0) {
//		printf( "Broadcasting dir\n");
		tree_broadcast_directory(dir);
	}
	return dir;
}

void tree::create_children() {
	auto boxl = tptr->box;
	auto boxr = tptr->box;
	const auto dim = tptr->level % NDIM;
	boxl.max[dim] = boxr.min[dim] = 0.5 * (tptr->box.min[dim] + tptr->box.max[dim]);
	const auto &localities = hpx_localities();
	const int min_level = bits_to_level(localities.size());
	const int child_level = tptr->level + 1;
	std::uint64_t il, ir;
	if (child_level <= min_level) {
		const auto total_nodes = (1 << child_level);
		const auto idl = (tptr->id << std::uint64_t(1));
		const auto idr = (tptr->id << std::uint64_t(1)) + std::uint64_t(1);
		il = (idl - total_nodes) * localities.size() / total_nodes;
		ir = (idr - total_nodes) * localities.size() / total_nodes;
	} else {
		il = hpx::get_locality_id();
		ir = hpx::get_locality_id();
	}
	auto futl = hpx::new_ < tree > (localities[il], (tptr->id << std::uint64_t(1)), boxl, tptr->level + 1);
	auto futr = hpx::new_ < tree > (localities[ir], (tptr->id << std::uint64_t(1)) + std::uint64_t(1), boxr, tptr->level + 1);
	tptr->leaf = false;
	auto id_left = futl.get();
	auto id_right = futr.get();
	hpx::future < std::uint64_t > fptrl;
	if (hpx::get_colocation_id(id_left).get() == hpx::find_here()) {
		fptrl = hpx::make_ready_future(get_ptr_action()(id_left));
	} else {
		fptrl = hpx::async < get_ptr_action > (id_left);
	}
	auto ptrr = get_ptr_action()(id_right);
	tptr->left_child = tree_client(std::move(id_left), fptrl.get());
	tptr->right_child = tree_client(std::move(id_right), ptrr);
}

int tree::destroy(int stack_cnt) {
	if (!tptr->leaf) {
		auto futl = tptr->left_child.destroy(stack_cnt);
		auto futr = tptr->right_child.destroy(stack_cnt);
		futl.get();
		tptr->left_child = tree_client();
		futr.get();
		tptr->right_child = tree_client();
	}
	if( tptr->level == 0 ) {
		tree_cleanup();
	}
	return 0;
}

check_pair tree::get_child_checks() const {
	check_pair checks;
	checks.first.info = &(tptr->child_info[0]);
	checks.first.opened = false;
	checks.second.info = &(tptr->child_info[1]);
	checks.second.opened = false;
	return std::move(checks);
}

std::uint64_t tree::get_ptr() {
	return reinterpret_cast<std::uint64_t>(this);
}

std::uint64_t tree::grow(int stack_cnt, bucket &&parts) {
	const int min_level = bits_to_level(hpx_localities().size());
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
		auto futl = tptr->left_child.grow(stack_cnt, true, std::move(parts_left));
		auto futr = tptr->right_child.grow(stack_cnt, false, std::move(parts_right));
		tptr->child_cnt[0] = futl.get();
		tptr->child_cnt[1] = futr.get();
	}
	return size();
}

bucket tree::get_parts() {
	return std::move(tptr->parts);
}

int tree::load_balance(int stack_cnt, std::uint64_t index, std::uint64_t total) {
	if (!tptr->leaf) {
		auto futl = tptr->left_child.load_balance(stack_cnt, true, index, total);
		auto futr = tptr->right_child.load_balance(stack_cnt, false, index + tptr->child_cnt[0], total);
		futl.get();
		futr.get();
		const auto &localities = hpx_localities();
		const int min_level = bits_to_level(localities.size());
		std::uint64_t il, ir;
		const int child_level = tptr->level + 1;
		if (child_level > min_level) {
			il = index * std::uint64_t(localities.size()) / total;
			ir = (index + tptr->child_cnt[0]) * std::uint64_t(localities.size()) / total;
			assert(il >= 0);
			assert(il < localities.size());
			assert(ir >= 0);
			assert(ir < localities.size());
			hpx::future<tree_client> newl;
			hpx::future<tree_client> newr;
			if (il != hpx::get_locality_id()) {
//				printf( "Migrating from %i to %i\n", hpx::get_locality_id(), il);
				newl = tptr->left_child.migrate(localities[il]);
			}
			if (ir != hpx::get_locality_id()) {
//				printf( "Migrating from %i to %i\n", hpx::get_locality_id(), ir);
				newr = tptr->right_child.migrate(localities[ir]);
			}
			if (newl.valid()) {
				tptr->left_child = newl.get();
			}
			if (newr.valid()) {
				tptr->right_child = newr.get();
			}
		}
	}
	return 0;
}

tree_client tree::migrate(hpx::id_type locality) {
	tptr->child_info[0].multi = nullptr;
	tptr->child_info[1].multi = nullptr;
	auto new_ptr = hpx::new_ < tree > (locality, *this);
	auto new_id = new_ptr.get();
	if (hpx::get_colocation_id(new_id).get() == hpx::find_here()) {
		return tree_client(new_id, get_ptr_action()(new_id));
	} else {
		return tree_client(new_id, hpx::async < get_ptr_action > (new_id).get());
	}

}

std::uint64_t tree::prune(int stack_cnt) {
	tptr->parent = tree_client();
	if (!tptr->leaf) {
		auto futl = tptr->left_child.prune(stack_cnt, true);
		auto futr = tptr->right_child.prune(stack_cnt, false);
		tptr->child_cnt[0] = futl.get();
		tptr->child_cnt[1] = futr.get();
		if (size() <= opts.bucket_size) {
			auto futl = tptr->left_child.get_parts();
			auto futr = tptr->right_child.get_parts();
			auto tmpl = futl.get();
			tptr->left_child = tree_client();
			auto tmpr = futr.get();
			tptr->right_child = tree_client();
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
	const int min_level = bits_to_level(hpx_localities().size());
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
		auto futl = tptr->left_child.verify(stack_cnt, true);
		auto futr = tptr->right_child.verify(stack_cnt, false);
		rc |= futl.get();
		rc |= futr.get();
	}
	return rc;
}

std::uint64_t tree::drift(int stack_cnt, int step, tree_client parent, tree_client self, float dt) {
	const int min_level = bits_to_level(hpx_localities().size());
	tptr->parent = parent;
	std::uint64_t drifted = 0;
	if (tptr->leaf) {
		while (tptr->lock++ != 0) {
			tptr->lock--;
		}
		bucket exit_parts;
		auto i = tptr->parts.begin();
		while (i != tptr->parts.end()) {
			if (step % 2 == i->step) {
				auto x = pos_to_double(i->x);
				const auto v = i->v;
				x += v * dt;
				i->x = double_to_pos(x);
				i->step++;
				if (!in_range(pos_to_double(i->x), tptr->box)) {
					exit_parts.insert(*i);
					drifted++;
					i = tptr->parts.remove(i);
				} else {
					i++;
				}
			} else {
				i++;
			}
		}
		tptr->lock--;
		if (exit_parts.size()) {
			if (tptr->level <= min_level) {
				directory->find_home(stack_cnt,std::move(exit_parts));
			} else {
				tptr->parent.find_home_parent(stack_cnt, std::move(exit_parts));
			}
		}
	} else {
		auto futl = tptr->left_child.drift(stack_cnt, true, step, self, tptr->left_child, dt);
		auto futr = tptr->right_child.drift(stack_cnt, false, step, self, tptr->right_child, dt);
		drifted += futl.get();
		drifted += futr.get();
	}
	if (tptr->level == 0) {
		tree_complete_drift();
	}
	return drifted;
}

int tree::find_home_child(int stack_cnt, bucket parts) {
	if (tptr->leaf) {
		while (tptr->lock++ != 0) {
			tptr->lock--;
		}
		while (parts.size()) {
			assert(in_range(pos_to_double(parts.front().x), tptr->box));
			tptr->parts.insert(parts.front());
			parts.remove(parts.begin());
		}
		tptr->lock--;
	} else {
		bucket l_parts;
		bucket r_parts;
		const auto dim = tptr->level % NDIM;
		const auto midx = 0.5 * (tptr->box.min[dim] + tptr->box.max[dim]);
		while (parts.size()) {
			auto &p = parts.front();
			const auto x = pos_to_double(p.x);
			if (x[dim] >= midx) {
				r_parts.insert(p);
			} else {
				l_parts.insert(p);
			}
			parts.remove(parts.begin());
		}
		if (l_parts.size()) {
			assert(tptr->left_child != tree_client());
			tptr->left_child.find_home_child(stack_cnt, std::move(l_parts));
		}
		if (r_parts.size()) {
			assert(tptr->right_child != tree_client());
			tptr->right_child.find_home_child(stack_cnt, std::move(r_parts));
		}
	}
	return 0;
}

int tree::find_home_parent(int stack_cnt, bucket parts) {
	const int min_level = bits_to_level(hpx_localities().size());
	bucket p_parts;
	bucket l_parts;
	bucket r_parts;
	const auto dim = tptr->level % NDIM;
	const auto midx = 0.5 * (tptr->box.min[dim] + tptr->box.max[dim]);
	while (parts.size()) {
		auto &p = parts.front();
		const auto x = pos_to_double(p.x);
		if (!in_range(x, tptr->box)) {
			p_parts.insert(p);
		} else if (x[dim] >= midx) {
			r_parts.insert(p);
		} else {
			l_parts.insert(p);
		}
		parts.remove(parts.begin());
	}
	if (p_parts.size()) {
		assert(tptr->parent != tree_client());
		if (tptr->level <= min_level) {
			directory->find_home(stack_cnt,std::move(p_parts));
		} else {
			tptr->parent.find_home_parent(stack_cnt, std::move(p_parts));
		}
	}
	if (l_parts.size()) {
		assert(tptr->left_child != tree_client());
		tptr->left_child.find_home_child(stack_cnt, std::move(l_parts));
	}
	if (r_parts.size()) {
		assert(tptr->right_child != tree_client());
		tptr->right_child.find_home_child(stack_cnt, std::move(r_parts));
	}
	return 0;
}
