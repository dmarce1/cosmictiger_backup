#include <cosmictiger/defs.hpp>
#include <cosmictiger/tree.hpp>

#include <atomic>
#include <stack>
#include <thread>

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<tree>, tree);

using destroy_action_type = tree::destroy_action;
using drift_action_type = tree::drift_action;
using drift_into_action_type = tree::drift_into_action;
using get_family_checks_action_type = tree::get_family_checks_action;
using find_family_action_type = tree::find_family_action;
using get_parts_action_type = tree::get_parts_action;
using get_ptr_action_type = tree::get_ptr_action;
using grow_action_type = tree::grow_action;
using prune_action_type = tree::prune_action;
using load_balance_action_type = tree::load_balance_action;
using migrate_action_type = tree::migrate_action;
using verify_action_type = tree::verify_action;

HPX_REGISTER_ACTION (destroy_action_type);
HPX_REGISTER_ACTION (drift_action_type);
HPX_REGISTER_ACTION (drift_into_action_type);
HPX_REGISTER_ACTION (get_family_checks_action_type);
HPX_REGISTER_ACTION (find_family_action_type);
HPX_REGISTER_ACTION (get_parts_action_type);
HPX_REGISTER_ACTION (get_ptr_action_type);
HPX_REGISTER_ACTION (grow_action_type);
HPX_REGISTER_ACTION (prune_action_type);
HPX_REGISTER_ACTION (load_balance_action_type);
HPX_REGISTER_ACTION (migrate_action_type);
HPX_REGISTER_ACTION (verify_action_type);

static std::stack<tree_mems*> stack;
static mutex_type mtx;

#define CHUNK_SIZE (512*1024)

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

tree_mems* tree_mems_alloc() {
	std::lock_guard < mutex_type > lock(mtx);
	if (stack.empty()) {
		tree_mems *base = (tree_mems*) malloc(sizeof(tree_mems) * CHUNK_SIZE);
		for (int i = 0; i < CHUNK_SIZE; i++) {
			stack.push(base + i);
		}
	}
	tree_mems *ptr = stack.top();
	stack.pop();
	return ptr;
}

void tree_mems_free(tree_mems *ptr) {
	std::lock_guard < mutex_type > lock(mtx);
	stack.push(ptr);
}

tree::tree() {
	tptr = tree_mems_alloc();
}

tree::tree(box_id_type id) {
	tptr = tree_mems_alloc();
	tptr->boxid = id;
	tptr->leaf = true;
}

tree::tree(const tree &other) {
	*tptr = *(other.tptr);
}

tree::~tree() {
	tree_mems_free(tptr);
}

void tree::create_children() {
	auto futl = hpx::new_ < tree > (hpx::find_here(), (tptr->boxid << 1) + 0);
	auto futr = hpx::new_ < tree > (hpx::find_here(), (tptr->boxid << 1) + 1);
	tptr->leaf = false;
	auto id_left = futl.get();
	auto id_right = futr.get();
	auto fptrl = hpx::async < get_ptr_action > (id_left);
	auto fptrr = hpx::async < get_ptr_action > (id_right);
	tptr->children[0] = tree_client(std::move(id_left), fptrl.get());
	tptr->children[1] = tree_client(std::move(id_right), fptrr.get());
}

void tree::destroy() {
	if (!tptr->leaf) {
		auto futl = tptr->children[0].destroy();
		auto futr = tptr->children[1].destroy();
		futl.get();
		futr.get();
	}
	tptr->self = tree_client();
}

int tree::drift(int stack_cnt, int step, float dt) {
	if (tptr->leaf) {
		std::array<bucket, NSIBLING> exit_parts;
		bucket parent_parts;
		{
			std::lock_guard < mutex_type > lock(tptr->mtx);
			const auto box = box_id_to_range(tptr->boxid);
			for (auto iter = tptr->parts.begin(); iter != tptr->parts.end();) {
				if (iter->step == step % 2) {
					auto &part = *iter;
					const auto v = part.v;
					const auto x = pos_to_double(part.x);
					int si = -1;
					part.dt = dt;
					auto this_dt = dt;
					for (int dim = dim; dim < NDIM; dim++) {
						if (v[dim] > 0.0) {
							const double t = (box.max[dim] - x[dim]) / v[dim];
							if (t > 0.0 && t < this_dt) {
								si = 2 * dim + 1;
								this_dt = t;
							}
						} else if (v[dim] < 0.0) {
							const double t = (box.min[dim] - x[dim]) / v[dim];
							if (t > 0.0 && t < this_dt) {
								si = 2 * dim;
								this_dt = t;
							}
						}
					}
					vect<double> this_x = x + v * this_dt;
					bool leaves_box = false;
					if (si != -1) {
						const auto dim = si / 2;
						part.dt -= this_dt;
						if (si % 2 == 0) {
							leaves_box = (this_x[dim] + v[dim] * part.dt) > box.max[dim];
						} else {
							leaves_box = (this_x[dim] + v[dim] * part.dt) < box.min[dim];
						}
						if (leaves_box) {
							part.x = this_x;
							if ((hpx::id_type) tptr->siblings[si] != hpx::invalid_id) {
								exit_parts[si].insert(part);
							} else {
								parent_parts.insert(part);
							}
							tptr->parts.remove(iter);
						}
					}
					if (!leaves_box) {
						part.x = this_x;
						part.dt = 0.0;
						part.step++;
						iter++;
					} else {
						printf("particle leaving box\n");
					}
				}
			}
		}
		std::array<hpx::future<void>, NSIBLING> futs;
		hpx::future<void> pfut;
		if (parent_parts.size()) {
			assert(tptr->parent != hpx::invalid_id);
			pfut = tptr->parent.drift_into(-1, std::move(parent_parts));
		}
		for (int si = 0; si < NSIBLING; si++) {
			if (exit_parts[si].size()) {
				futs[si] = tptr->siblings[si].drift_into(si ^ 1, std::move(exit_parts[si]));
			} else {
				futs[si] = hpx::make_ready_future();
			}
		}
		if (pfut.valid()) {
			pfut.get();
		}
		hpx::wait_all(futs.begin(), futs.end());
	} else {
		auto futl = tptr->children[0].drift(stack_cnt, step, dt);
		auto futr = tptr->children[1].drift(stack_cnt, step, dt);
		futl.get();
		futr.get();
	}
	return 0;
}

void tree::drift_into(int sender, bucket parts) {
	std::array<bucket, NSIBLING> exit_parts;
	bucket parent_parts;
	{
		std::lock_guard < mutex_type > lock(tptr->mtx);
		const auto box = box_id_to_range(tptr->boxid);
		for (auto iter = parts.begin(); iter != parts.end();) {
			auto &part = *iter;
			const auto v = part.v;
			const auto x = pos_to_double(part.x);
			int si = -1;
			auto this_dt = part.dt;
			for (int dim = dim; dim < NDIM; dim++) {
				if (v[dim] > 0.0) {
					const double t = (box.max[dim] - x[dim]) / v[dim];
					if (t > 0.0 && t < this_dt) {
						si = 2 * dim + 1;
						this_dt = t;
					}
				} else if (v[dim] < 0.0) {
					const double t = (box.min[dim] - x[dim]) / v[dim];
					if (t > 0.0 && t < this_dt) {
						si = 2 * dim;
						this_dt = t;
					}
				}
			}
			vect<double> this_x = x + v * this_dt;
			bool leaves_box = false;
			if (si != -1) {
				const auto dim = si / 2;
				part.dt -= this_dt;
				if (si % 2 == 0) {
					leaves_box = (this_x[dim] + v[dim] * part.dt) > box.max[dim];
				} else {
					leaves_box = (this_x[dim] + v[dim] * part.dt) < box.min[dim];
				}
				if (leaves_box) {
					part.x = this_x;
					if ((hpx::id_type) tptr->siblings[si] != hpx::invalid_id) {
						exit_parts[si].insert(part);
					} else {
						parent_parts.insert(part);
					}
					tptr->parts.remove(iter);
				}
			}
			if (!leaves_box) {
				part.x = this_x;
				part.dt = 0.0;
				part.step++;
				tptr->parts.insert(part);
				iter++;
			}
		}
	}
	hpx::future < std::uint64_t > grow_fut;
	if (tptr->parts.size() > opts.bucket_size) {
		grow_fut = hpx::async([this](bucket &&parts) {
			return grow(0, std::move(parts));
		}, std::move(tptr->parts));
	}
	std::array<hpx::future<void>, NSIBLING> futs;
	hpx::future<void> pfut;
	if (parent_parts.size()) {
		assert(tptr->parent != hpx::invalid_id);
		pfut = tptr->parent.drift_into(-1, std::move(parent_parts));
	}
	for (int si = 0; si < NSIBLING; si++) {
		if (exit_parts[si].size()) {
			futs[si] = tptr->siblings[si].drift_into(si ^ 1, std::move(exit_parts[si]));
		} else {
			futs[si] = hpx::make_ready_future();
		}
	}
	if (pfut.valid()) {
		pfut.get();
	}
	if (grow_fut.valid()) {
		grow_fut.get();
	}
	hpx::wait_all(futs.begin(), futs.end());
}

std::array<family_check, NCHILD> tree::get_family_checks() const {
	std::array<family_check, NCHILD> checks;
	for (int ci = 0; ci < NCHILD; ci++) {
		checks[ci].node = tptr->children[ci];
		checks[ci].boxid = (tptr->boxid << 1) + ci;
	}
	return checks;
}

std::uint64_t tree::get_ptr() {
	return reinterpret_cast<std::uint64_t>(this);
}

int tree::find_family(int stack_cnt, tree_client parent, tree_client self, std::vector<family_check> checks) {
	tptr->parent = parent;
	tptr->self = self;
	for (int si = 0; si < NSIBLING; si++) {
		tptr->siblings[si] = tree_client();
	}
	const auto box = box_id_to_range(tptr->boxid);
	for (const auto &c : checks) {
		const auto obox = box_id_to_range(c.boxid);
		const auto si = range_sibling_index(box, obox);
		if (si != -1) {
			tptr->siblings[si] = c.node;
		}
	}
	if (!tptr->leaf) {
		std::vector < hpx::future<std::array<family_check, NCHILD>> > cfuts;
		for (const auto &c : checks) {
			cfuts.push_back(c.node.get_family_checks());
		}
		checks.resize(0);
		checks.reserve(cfuts.size() * NCHILD);
		for (auto &fut : cfuts) {
			const auto tmp = fut.get();
			for (int ci = 0; ci < NCHILD; ci++) {
				checks.push_back(tmp[ci]);
			}
		}
		auto futl = tptr->children[0].find_family(stack_cnt, self, tptr->children[0], checks);
		auto futr = tptr->children[1].find_family(stack_cnt, self, tptr->children[1], std::move(checks));
		futl.get();
		futr.get();
	}
	return 0;
}

std::uint64_t tree::grow(int stack_cnt, bucket parts) {

	tptr->child_cnt[0] = 0;
	tptr->child_cnt[1] = 0;
	if (tptr->leaf) {
		if (parts.size() + tptr->parts.size() <= opts.bucket_size) {
			while (parts.size()) {
				tptr->parts.insert(parts.front());
				parts.remove(parts.begin());
			}
		} else {
			create_children();
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
		auto futl = tptr->children[0].grow(stack_cnt, std::move(parts_left));
		auto futr = tptr->children[1].grow(stack_cnt, std::move(parts_right));
		tptr->child_cnt[0] = futl.get();
		tptr->child_cnt[1] = futr.get();
	}
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
		const std::uint64_t il = index * opts.problem_size / localities.size();
		const std::uint64_t ir = (index + tptr->child_cnt[0]) * opts.problem_size / localities.size();
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
		auto futl = tptr->children[0].load_balance(stack_cnt, il);
		auto futr = tptr->children[1].load_balance(stack_cnt, ir);
		futl.get();
		futr.get();
	}
	return 0;
}

tree_client tree::migrate(hpx::id_type locality) {
	auto new_ptr = hpx::new_ < tree > (locality, *this);
	tptr->self = tree_client(new_ptr.get(), this);
	return tptr->self;

}

int tree::prune(int stack_cnt) {
	if (!tptr->leaf) {
		if (size() <= opts.bucket_size) {
			auto futl = tptr->children[0].get_parts();
			auto futr = tptr->children[1].get_parts();
			auto tmp = futl.get();
			while (tmp.size()) {
				tptr->parts.insert(tmp.front());
				tmp.remove(tmp.begin());
			}
			tmp = futr.get();
			while (tmp.size()) {
				tptr->parts.insert(tmp.front());
				tmp.remove(tmp.begin());
			}
			tptr->leaf = true;
		}
	} else {
		auto futl = tptr->children[0].prune(stack_cnt);
		auto futr = tptr->children[1].prune(stack_cnt);
		futl.get();
		futr.get();
	}
	return 0;
}

std::size_t tree::size() const {
	if (tptr->leaf) {
		return tptr->parts.size();
	} else {
		return tptr->child_cnt[0] + tptr->child_cnt[1];
	}
}

int tree::verify(int stack_cnt) const {
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
		if (size() <= opts.bucket_size) {
			printf("%i %i\n", tptr->child_cnt[0], tptr->child_cnt[1]);
			rc |= TREE_UNDERFLOW;
		}
		auto futl = tptr->children[0].verify(stack_cnt);
		auto futr = tptr->children[1].verify(stack_cnt);
		rc |= futl.get();
		rc |= futr.get();
	}
	return rc;
}

