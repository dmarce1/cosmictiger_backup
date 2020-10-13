#include <cosmictiger/defs.hpp>
#include <cosmictiger/util.hpp>
#include <cosmictiger/tree.hpp>
#include <cosmictiger/tree_dir.hpp>
#include <cosmictiger/gravity.hpp>
#include <cosmictiger/gravity_queue.hpp>
#include <cosmictiger/pos_cache.hpp>
#include <cosmictiger/output.hpp>
#include <cosmictiger/time.hpp>

#include <atomic>
#include <stack>
#include <thread>

#define WORKGROUP_SIZE 64

#ifdef HPX_LITE
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<tree>, tree);
using build_tree_dir_type = tree::build_tree_dir_action;
using compute_multipoles_type = tree::compute_multipoles_action;
using destroy_type = tree::destroy_action;
using drift_type = tree::drift_action;
using place_parts_type = tree::place_parts_action;
using grow_type = tree::grow_action;
using kick_fmm_type = tree::kick_fmm_action;
using load_balance_type = tree::load_balance_action;
using count_children_type = tree::count_children_action;
using verify_type = tree::verify_action;
using get_child_checks_type = tree::get_child_checks_action;
using get_check_info_type = tree::get_check_info_action;
using get_ptr_type = tree::get_ptr_action;
using get_parts_type = tree::get_parts_action;
using get_positions_type = tree::get_positions_action;
using migrate_type = tree::migrate_action;
HPX_REGISTER_ACTION(build_tree_dir_type);
HPX_REGISTER_ACTION(compute_multipoles_type);
HPX_REGISTER_ACTION(destroy_type);
HPX_REGISTER_ACTION(drift_type);
HPX_REGISTER_ACTION(place_parts_type);
HPX_REGISTER_ACTION(grow_type);
HPX_REGISTER_ACTION(kick_fmm_type);
HPX_REGISTER_ACTION(load_balance_type);
HPX_REGISTER_ACTION(count_children_type);
HPX_REGISTER_ACTION(verify_type);
HPX_REGISTER_ACTION(get_child_checks_type);
HPX_REGISTER_ACTION(get_check_info_type);
HPX_REGISTER_ACTION(get_ptr_type);
HPX_REGISTER_ACTION(get_parts_type);
HPX_REGISTER_ACTION(get_positions_type);
HPX_REGISTER_ACTION(migrate_type);

#else
HPX_REGISTER_COMPONENT(hpx::components::managed_component<tree>, tree);
#endif

static std::shared_ptr<tree_dir> directory;
static fmm_params fmm;
static kick_return krc;
static mutex_type krc_mtx;

void tree_broadcast_directory(const tree_dir &dir);
void tree_cleanup();

HPX_PLAIN_ACTION (tree_broadcast_directory);
HPX_PLAIN_ACTION (tree_cleanup);
HPX_PLAIN_ACTION (tree_set_fmm_params);

tree_client tree_allocate(tree node) {
	auto id = hpx::new_<tree>(hpx::find_here(), std::move(node)).get();
	auto ptr = tree::get_ptr_action()(id);
	return tree_client(id, ptr);
}

HPX_PLAIN_ACTION (tree_allocate);

int tree_ewald_min_level(double theta, double h) {
	if (opts.ewald) {
		int lev = 12;
		while (1) {
			int N = 1 << (lev / NDIM);
			double dx = 0.25 * N;
			double a;
			if (lev % NDIM == 0) {
				a = std::sqrt(3);
			} else if (lev % NDIM == 1) {
				a = 1.5;
			} else {
				a = std::sqrt(1.5);
			}
			double r = (1.0 + opts.sink_bias) * a / theta + h * N;
			if (dx > r) {
				break;
			}
			lev++;
		}
		static std::atomic<int> lk(0);
		if (lk++ == 0 && hpx::get_locality_id() == 0) {
			printf("Ewald level = %i\n", lev);
		}
		return lev;
	} else {
		return 0;
	}
}

void tree_insert_parts(bucket &&parts) {
	directory->find_home(std::move(parts)).get();
}

kick_return tree_kick_return();

HPX_PLAIN_ACTION (tree_kick_return);

kick_return tree_kick_return() {
	const auto &localities = hpx_localities();
	const int il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const int ir = ((hpx::get_locality_id() + 1) << 1);
	std::vector < hpx::future < kick_return >> futs;
	if (il < localities.size()) {
		futs.push_back(hpx::async<tree_kick_return_action>(localities[il]));
	}
	if (ir < localities.size()) {
		futs.push_back(hpx::async<tree_kick_return_action>(localities[ir]));
	}

	for (int i = 0; i < futs.size(); i++) {
		auto tmp = futs[i].get();
		krc.rung = std::max(krc.rung, tmp.rung);
		krc.epot += tmp.epot;
	}
	kick_return rc = krc;
	krc = kick_return();
	return rc;
}

void tree_broadcast_directory(const tree_dir &dir) {
	const auto &localities = hpx_localities();
	const int il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const int ir = ((hpx::get_locality_id() + 1) << 1);
	std::vector<hpx::future<void>> futs;
	if (il < localities.size()) {
		futs.push_back(hpx::async<tree_broadcast_directory_action>(localities[il], dir));
	}
	if (ir < localities.size()) {
		futs.push_back(hpx::async<tree_broadcast_directory_action>(localities[ir], dir));
	}
//	printf( "Setting dir %i\n", dir.size());
	directory = std::make_shared < tree_dir > (dir);
	hpx::wait_all(futs.begin(), futs.end());
}

void tree_set_fmm_params(fmm_params p) {
	const auto &localities = hpx_localities();
	const int il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const int ir = ((hpx::get_locality_id() + 1) << 1);
	std::vector<hpx::future<void>> futs;
	if (il < localities.size()) {
		futs.push_back(hpx::async<tree_set_fmm_params_action>(localities[il], p));
	}
	if (ir < localities.size()) {
		futs.push_back(hpx::async<tree_set_fmm_params_action>(localities[ir], p));
	}
	fmm = p;
	hpx::wait_all(futs.begin(), futs.end());
}

void tree_cleanup() {
	const auto &localities = hpx_localities();
	const int il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const int ir = ((hpx::get_locality_id() + 1) << 1);
	std::vector<hpx::future<void>> futs;
	if (il < localities.size()) {
		futs.push_back(hpx::async<tree_cleanup_action>(localities[il]));
	}
	if (ir < localities.size()) {
		futs.push_back(hpx::async<tree_cleanup_action>(localities[ir]));
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
	tptr->child_cnt[0] = 0;
	tptr->child_cnt[1] = 0;
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
		auto futl = hpx::async<build_tree_dir_action>(tptr->children[0].get_id(), tptr->children[0]);
		auto futr = hpx::async<build_tree_dir_action>(tptr->children[1].get_id(), tptr->children[1]);
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

multipole_return tree::compute_multipoles(int stack_cnt, std::uint64_t work_id, bucket &&parts, std::uint64_t index) {

	const auto dist_level = bits_to_level(hpx_localities().size());
	const int min_level = std::max(dist_level, tree_ewald_min_level(fmm.theta, opts.h));
	multipole_return rc;
	std::uint64_t nactive = 0;
	vect<double> xc;
	multipole<float> M;
	M = float(0.0);
	range prange;
	float r = 0.0;
	for (int dim = 0; dim < NDIM; dim++) {
		xc[dim] = 0.0;
		prange.max[dim] = 0.0;
		prange.min[dim] = 1.0;
	}
	if (work_id == -1 && size() <= WORKGROUP_SIZE * opts.bucket_size) {
		work_id = gravity_queue_genid();
	}
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
	} else if (tptr->level > min_level) {
		const auto sz = parts.size() + tptr->parts.size() + tptr->child_cnt[0] + tptr->child_cnt[1];
		if (sz <= opts.bucket_size) {
			for (auto i = parts.begin(); i != parts.end(); i++) {
				tptr->parts.insert(*i);
			}
			auto tmp = get_parts();

			tptr->parts = std::move(tmp);
			tptr->children[0] = tree_client();
			tptr->children[1] = tree_client();
			tptr->child_cnt[0] = 0;
			tptr->child_cnt[1] = 0;
			tptr->leaf = true;
		}
	}
	if (tptr->leaf) {
		const auto m = opts.particle_mass;
		if (tptr->parts.size()) {
			for (auto i = tptr->parts.begin(); i != tptr->parts.end(); i++) {
				const auto x = pos_to_double(i->x);
				if (fmm.stats || (i->rung >= fmm.min_rung)) {
					nactive++;
				}
				for (int dim = 0; dim < NDIM; dim++) {
					xc[dim] += x[dim];
					prange.max[dim] = std::max(prange.max[dim], x[dim]);
					prange.min[dim] = std::min(prange.min[dim], x[dim]);
				}
			}
			xc /= tptr->parts.size();
			for (auto i = tptr->parts.begin(); i != tptr->parts.end(); i++) {
				const auto x = pos_to_double(i->x);
				r = std::max(r, (float) abs(x - xc));
				M() += m;
				for (int n = 0; n < NDIM; n++) {
					const auto dxn = x[n] - xc[n];
					for (int p = 0; p <= n; p++) {
						const auto dxp = x[p] - xc[p];
						M(n, p) += dxn * dxp * m;
						for (int l = 0; l <= p; l++) {
							const auto dxl = x[l] - xc[l];
							M(n, p, l) -= dxn * dxp * dxl * m;
						}
					}
				}
			}
		} else {
			xc = range_center(tptr->box);
			prange.max = prange.min = xc;
		}
	} else {
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
		tptr->child_cnt[0] += parts_left.size();
		tptr->child_cnt[1] += parts_right.size();
		local_load_balance(index, opts.problem_size);
		auto futl = tptr->children[0].compute_multipoles(stack_cnt, true, work_id, std::move(parts_left), index);
		auto futr = tptr->children[1].compute_multipoles(stack_cnt, false, work_id, std::move(parts_right), index + tptr->child_cnt[0]);
		auto L = futl.get();
		auto R = futr.get();
		tptr->child_info[0] = L.info;
		tptr->child_info[1] = R.info;
		nactive = L.nactive + R.nactive;
		const multi_src &ml = *L.info.multi;
		const multi_src &mr = *R.info.multi;
		const auto mtot = ml.m() + mr.m();
		if (mtot != 0.0) {
			const auto xl = pos_to_double(ml.x);
			const auto xr = pos_to_double(mr.x);
			xc = (xl * ml.m() + xr * mr.m()) / mtot;
			const auto dxl = xl - xc;
			const auto dxr = xr - xc;
			M = (ml.m >> dxl) + (mr.m >> dxr);
			if (mr.m() == 0.0) {
				r = ml.r;
				prange = L.prange;
			} else if (ml.m() == 0.0) {
				r = mr.r;
				prange = R.prange;
			} else {
				r = std::max(abs(dxl) + ml.r, abs(dxr) + mr.r);
				for (int dim = 0; dim < NDIM; dim++) {
					prange.max[dim] = std::max(L.prange.max[dim], R.prange.max[dim]);
					prange.min[dim] = std::min(L.prange.min[dim], R.prange.min[dim]);
				}
			}
			vect<double> corner;
			float rmax = 0.0;
			for (int ci = 0; ci < 8; ci++) {
				for (int dim = 0; dim < NDIM; dim++) {
					corner[dim] = ((ci >> dim) & 1) ? prange.min[dim] : prange.max[dim];
				}
				rmax = std::max(rmax, (float) abs(corner - xc));
			}
			r = std::min(rmax, r);
		} else {
			xc = range_center(tptr->box);
			prange.max = prange.min = xc;
		}
	}
	tptr->multi.m = M;
	tptr->multi.x = double_to_pos(xc);
	assert(in_range(pos_to_double(tptr->multi.x), tptr->box));
	tptr->multi.r = r;
	tptr->nactive = nactive;
	tptr->work_id = work_id;
	if (tptr->leaf && tptr->nactive) {
		assert(work_id != -1);
		gravity_queue_checkin(work_id);
	}
	rc.info.multi = &tptr->multi;
	rc.info.leaf = tptr->leaf;
	rc.info.node.rank = hpx::get_locality_id();
	rc.info.node.ptr = reinterpret_cast<std::uint64_t>(this);
	rc.nactive = nactive;
	rc.prange = prange;
	return rc;
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
	auto futl = hpx::new_<tree>(localities[il], (tptr->id << std::uint64_t(1)), boxl, tptr->level + 1);
	auto futr = hpx::new_<tree>(localities[ir], (tptr->id << std::uint64_t(1)) + std::uint64_t(1), boxr, tptr->level + 1);
	tptr->leaf = false;
	auto id_left = futl.get();
	auto id_right = futr.get();
	hpx::future < std::uint64_t > fptrl;
	if (hpx::get_colocation_id(id_left).get() == hpx::find_here()) {
		fptrl = hpx::make_ready_future(get_ptr_action()(id_left));
	} else {
		fptrl = hpx::async<get_ptr_action>(id_left);
	}
	auto ptrr = get_ptr_action()(id_right);
	tptr->children[0] = tree_client(std::move(id_left), fptrl.get());
	tptr->children[1] = tree_client(std::move(id_right), ptrr);
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
	if (tptr->level == 0) {
		tree_cleanup();
	}
	return 0;
}

check_pair tree::get_child_checks() const {
	assert(!tptr->leaf);
	check_pair checks;
	checks.first.info = &(tptr->child_info[0]);
	checks.first.opened = false;
	checks.second.info = &(tptr->child_info[1]);
	checks.second.opened = false;
	return std::move(checks);
}

std::vector<part_pos> tree::get_positions() const {
	std::vector<part_pos> pos;
	pos.reserve(tptr->parts.size());
	for (auto i = tptr->parts.begin(); i != tptr->parts.end(); i++) {
		pos.push_back(i->x);
	}
	return std::move(pos);
}

check_info tree::get_check_info() const {
	check_info info;
	info.leaf = tptr->leaf;
	info.multi = &tptr->multi;
	info.node.rank = hpx::get_locality_id();
	info.node.ptr = reinterpret_cast<std::uint64_t>(this);
	return info;

}

std::uint64_t tree::get_ptr() {
	return reinterpret_cast<std::uint64_t>(this);
}

std::uint64_t tree::grow(int stack_cnt, bucket &&parts, bool first_pass) {
	const auto dist_level = bits_to_level(hpx_localities().size());
	const int min_level = first_pass ? dist_level : std::max(dist_level, tree_ewald_min_level(fmm.theta, opts.h));
//	printf("%i %i %i %i\n", tptr->level, min_level, (parts.size() + tptr->parts.size()), opts.bucket_size);
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
			//		printf("%li %e %e\n", i->x[dim], double(i->x[dim]), xmid);
			if (double(i->x[dim]) >= xmid) {
				parts_right.insert(*i);
				i = parts_left.remove(i);
			} else {
				i++;
			}
		}
		auto futl = tptr->children[0].grow(stack_cnt, true, std::move(parts_left), first_pass);
		auto futr = tptr->children[1].grow(stack_cnt, false, std::move(parts_right), first_pass);
		tptr->child_cnt[0] = futl.get();
		tptr->child_cnt[1] = futr.get();
	} else {
		tptr->child_cnt[0] = 0;
		tptr->child_cnt[1] = 0;
	}
//	if( tptr->child_cnt[0] > opts.problem_size) {
//		printf( "--------%li %i\n",tptr->child_cnt[1], tptr->leaf);
//	}
	return size();
}

bucket tree::get_parts() {
	if (!tptr->leaf) {
		auto futl = tptr->children[0].get_parts();
		auto futr = tptr->children[1].get_parts();
		auto L = futl.get();
		auto R = futr.get();
		while (L.size()) {
			tptr->parts.insert(L.front());
			L.remove(L.begin());
		}
		while (R.size()) {
			tptr->parts.insert(R.front());
			R.remove(R.begin());
		}
	}
	return std::move(tptr->parts);
}

std::vector<int> tree::checks_far(const std::vector<check_item> &checks, bool ewald) {
	static const float h = opts.h;
	const simd_float R1 = tptr->multi.r * opts.sink_bias + h * fmm.theta;
	vect<simd_int> X1;
	for (int dim = 0; dim < NDIM; dim++) {
		X1[dim] = (int) tptr->multi.x[dim];
	}
	const simd_float Thetainv = 1.0 / fmm.theta;
	const int simd_size = (checks.size() - 1) / simd_float::size() + 1;
	std::vector<simd_float> R2(simd_size);
	std::vector<vect<simd_int>> X2(simd_size);
	std::vector<int> res(checks.size());
	for (int i = 0; i < checks.size(); i += simd_float::size()) {
		const int j = i / simd_float::size();
		for (int k = 0; k < simd_float::size(); k++) {
			const int l = std::min(i + k, (int) checks.size() - 1);
			R2[j][k] = checks[l].info->multi->r;
			for (int dim = 0; dim < NDIM; dim++) {
				X2[j][dim][k] = (int) checks[l].info->multi->x[dim];
			}
		}
	}
	for (int i = 0; i < simd_size; i++) {
		vect<simd_float> dX;
		if (opts.ewald) {
			const vect<simd_int> dXi = X2[i] - X1;
			for (int dim = 0; dim < NDIM; dim++) {
				dX[dim] = simd_float(dXi[dim]) * POS_INV;
			}
		} else {
			for (int dim = 0; dim < NDIM; dim++) {
				dX[dim] = (simd_float(X2[i][dim]) - simd_float(X1[dim])) * POS_INV;
			}
		}
		simd_float dist = abs(dX);
		if (ewald) {
			dist = max(simd_float(0.25), dist);
		}
		const simd_float b = dist > (R2[i] + R1) * Thetainv;
		for (int j = 0; j < simd_float::size(); j++) {
			const auto k = std::min(simd_float::size() * i + j, res.size() - 1);
			res[k] = b[j];
		}
	}
	return res;
}

int tree::kick_fmm(int stack_cnt, std::vector<check_item> &&dchecks, std::vector<check_item> &&echecks, expansion_src &&L) {
#define LOAD_SIZE (16 * simd_float::size())
	std::vector<check_item> next_dchecks;
	std::vector<check_item> next_echecks;
	next_dchecks.reserve(dchecks.size());
	next_echecks.reserve(echecks.size());
	std::vector<tree_ptr> PP_list;
	std::vector<tree_ptr> CP_list;
	std::vector<const multi_src*> PC_list;
	std::vector<const multi_src*> CC_list;
	std::vector<const multi_src*> ewald_list;
	const auto xcom = pos_to_double(tptr->multi.x);
	if (tptr->nactive > 0) {
		L.l = L.l << (xcom - L.x);
		L.x = xcom;
		auto far = checks_far(dchecks, false);
		for (int i = 0; i < dchecks.size(); i++) {
			const auto &pos = dchecks[i].info->node;
			const auto &mpole = dchecks[i].info->multi;
			if (far[i]) {
				if (dchecks[i].opened) {
					CP_list.push_back(pos);
				} else {
					CC_list.push_back(mpole);
					if (CC_list.size() == LOAD_SIZE) {
						gravity_CC_direct(L.l, tptr->multi.x, CC_list, fmm.stats);
						CC_list.resize(0);
					}
				}
			} else {
				next_dchecks.push_back(dchecks[i]);
			}
		}
		gravity_CC_direct(L.l, tptr->multi.x, CC_list, fmm.stats);
		if (opts.ewald) {
			far = checks_far(echecks, true);
			for (int i = 0; i < echecks.size(); i++) {
				const auto &pos = echecks[i].info->node;
				const auto &mpole = echecks[i].info->multi;
				if (far[i]) {
					assert(!echecks[i].opened);
					ewald_list.push_back(mpole);
				} else {
					next_echecks.push_back(echecks[i]);
					if (next_echecks.size() == LOAD_SIZE) {
						gravity_CC_ewald(L.l, tptr->multi.x, ewald_list, fmm.stats);
						next_echecks.resize(0);
					}

				}
			}
			gravity_CC_ewald(L.l, tptr->multi.x, ewald_list, fmm.stats);
		}

		hpx::future<std::vector<check_item>> echecks_fut, dchecks_fut;
		dchecks_fut = get_next_checklist(next_dchecks);
		if (opts.ewald) {
			echecks_fut = get_next_checklist(next_echecks);
		}
		auto cp_ptrs = ::get_positions(CP_list);
		static thread_local std::vector<part_pos> cp;
		cp.resize(0);
		for (int i = 0; i < cp_ptrs.size(); i++) {
			const auto &this_x = *cp_ptrs[i];
			for (int j = 0; j < this_x.size(); j++) {
				cp.push_back(this_x[j]);
				if (cp.size() == LOAD_SIZE) {
					gravity_CP_direct(L.l, tptr->multi.x, cp, fmm.stats);
					cp.resize(0);
				}
			}
		}
		gravity_CP_direct(L.l, tptr->multi.x, cp, fmm.stats);

		dchecks = dchecks_fut.get();
		if (opts.ewald) {
			echecks = echecks_fut.get();
		}
		if (tptr->leaf) {
			assert(!opts.ewald || echecks.size() == 0);
			while (dchecks.size()) {
				next_dchecks.resize(0);
				const auto far = checks_far(dchecks, false);
				for (int i = 0; i < dchecks.size(); i++) {
					const auto &pos = dchecks[i].info->node;
					const auto &mpole = dchecks[i].info->multi;
					if (dchecks[i].opened) {
						PP_list.push_back(pos);
					} else {
						if (far[i]) {
							PC_list.push_back(mpole);
						} else {
							next_dchecks.push_back(dchecks[i]);
						}
					}
				}
				dchecks_fut = get_next_checklist(next_dchecks);
				dchecks = dchecks_fut.get();
			}
			auto x = std::make_shared<std::vector<part_pos>>();
			for (auto i = tptr->parts.begin(); i != tptr->parts.end(); i++) {
				if (i->rung >= fmm.min_rung || fmm.stats) {
					x->push_back(i->x);
				}
			}
			auto f = std::make_shared < std::vector < _4force >> (x->size());
			for (auto i = 0; i != x->size(); i++) {
				_4force this_f;
				L.l.translate_L2(this_f.g, this_f.phi, vect<float>(vect<double>(pos_to_double((*x)[i]) - xcom)));
				(*f)[i].phi = this_f.phi;
				(*f)[i].g = this_f.g;
			}
			gravity_queue_add_work(tptr->work_id, f, x, std::move(PP_list), std::move(PC_list), [f, x, this]() {
				for (auto i = f->begin(); i != f->end(); i++) {
					i->phi += SELF_PHI * opts.particle_mass / opts.h;
					i->g = i->g * opts.G;
					i->phi = i->phi * opts.G;
				}
				int j = 0;
				int max_rung = 0;
				const auto ainv = 1.0 / fmm.a;
				const auto a3inv = ainv * ainv * ainv;
				const auto m = opts.particle_mass;
				double epot = 0.0;
				for (auto i = tptr->parts.begin(); i != tptr->parts.end(); i++) {
					if (fmm.stats || i->rung >= fmm.min_rung) {
						if (fmm.stats) {
							if (i->out) {
								output_add_particle(*i, (*f)[j]);
							}
							epot += 0.5 * m * (*f)[j].phi * ainv;
						}
						double this_dt;
						if (i->rung >= fmm.min_rung) {
							if (i->rung != 0) {
								this_dt = rung_to_dt(i->rung);
								i->v = i->v + (*f)[j].g * this_dt * ainv * 0.5;
							}
							const auto g = abs((*f)[j].g) * a3inv;
							this_dt = opts.eta * std::sqrt(opts.h / (g * SELF_PHI));
							i->rung = std::max(std::max((int) dt_to_rung(this_dt), (int) fmm.min_rung), (int) (i->rung - 1));
							max_rung = std::max((int) max_rung, (int) i->rung);
							this_dt = rung_to_dt(i->rung);
							i->v = i->v + (*f)[j].g * this_dt * ainv * 0.5;
						}
						j++;
					}
				}
				std::lock_guard<mutex_type> lock(krc_mtx);
				krc.rung = std::max(krc.rung, max_rung);
				if (fmm.stats) {
					krc.epot += epot;
				}
			}, fmm.stats);

		} else {
			auto futl = tptr->children[0].kick_fmm(stack_cnt, true, dchecks, echecks, L);
			auto futr = tptr->children[1].kick_fmm(stack_cnt, false, std::move(dchecks), std::move(echecks), L);
			futl.get();
			futr.get();
		}

	}
	if (tptr->level == 0) {
		check_cleanup();
		pos_cache_cleanup();
	}
	return 0;
}

std::uint64_t tree::local_load_balance(std::uint64_t index, std::uint64_t total) {
	std::uint64_t nmig = 0;
	const auto &localities = hpx_localities();
	const int min_level = bits_to_level(localities.size());
	std::uint64_t il, ir;
	const int child_level = tptr->level + 1;
	if (child_level > min_level) {
		il = index * std::uint64_t(localities.size()) / (total + 1);
		ir = (index + tptr->child_cnt[0]) * std::uint64_t(localities.size()) / (total + 1);
		assert(il >= 0);
		assert(il < localities.size());
		assert(ir >= 0);
		assert(ir < localities.size());
		hpx::future<tree_client> newl;
		hpx::future<tree_client> newr;
		if (il != tptr->children[0].get_rank()) {
			nmig++;
			newl = tptr->children[0].migrate(localities[il]);
		}
		if (ir != tptr->children[1].get_rank()) {
			nmig++;
			newr = tptr->children[1].migrate(localities[ir]);
		}
		if (newl.valid()) {
			tptr->children[0] = newl.get();
		}
		if (newr.valid()) {
			tptr->children[1] = newr.get();
		}
	}
	return nmig;
}

tree_stats tree::load_balance(int stack_cnt, std::uint64_t index, std::uint64_t total) {
	tree_stats stats;
	assert(index <= total);
	if (!tptr->leaf) {
		stats.nnode++;
		stats.nmig = local_load_balance(index, total);
		auto futl = tptr->children[0].load_balance(stack_cnt, true, index, total);
		auto futr = tptr->children[1].load_balance(stack_cnt, false, index + tptr->child_cnt[0], total);
		stats += futl.get();
		stats += futr.get();
	} else {
		stats.nnode++;
		stats.nleaf++;
	}
	return stats;
}

tree_client tree::migrate(hpx::id_type locality) {
	tptr->child_info[0].multi = nullptr;
	tptr->child_info[1].multi = nullptr;
	return tree_allocate_action()(locality, std::move(*this));
}

std::uint64_t tree::count_children() {
//	printf("count_children\n");
	const int min_level = bits_to_level(hpx_localities().size());
	if (tptr->level < min_level) {
		auto futl = tptr->children[0].count_children();
		auto futr = tptr->children[1].count_children();
		tptr->child_cnt[0] = futl.get();
		tptr->child_cnt[1] = futr.get();
	}
	return tptr->parts.size() + tptr->child_cnt[0] + tptr->child_cnt[1];
}

std::size_t tree::size() const {
	if (tptr->leaf) {
		return tptr->parts.size();
	} else {
		return tptr->child_cnt[0] + tptr->child_cnt[1];
	}
}

int tree::verify(int stack_cnt) const {
	const int min_level = std::max(bits_to_level(hpx_localities().size()), tree_ewald_min_level(fmm.theta, opts.h));
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
		auto futl = tptr->children[0].verify(stack_cnt, true);
		auto futr = tptr->children[1].verify(stack_cnt, false);
		rc |= futl.get();
		rc |= futr.get();
	}
	return rc;
}

drift_return tree::drift(int stack_cnt, double dt, double abar) {
	const auto dist_level = bits_to_level(hpx_localities().size());
	drift_return rc;
	bucket &exit_parts = rc.parts;

	rc.ekin = 0.0;
	const auto m = opts.particle_mass;
	if (tptr->leaf) {
		std::lock_guard<mutex_type> lock(tptr->mtx);
		auto i = tptr->parts.begin();
		while (i != tptr->parts.end()) {
			auto x = pos_to_double(i->x);
			const auto v = i->v;
			x += v * dt / (abar * abar);
			i->x = double_to_pos(x);
			rc.ekin += 0.5 * m * (v.dot(v) / (abar * abar));
			if (!in_range(pos_to_double(i->x), tptr->box)) {
				exit_parts.insert(*i);
				i = tptr->parts.remove(i);
			} else {
				i++;
			}
		}
		rc.ndrift = exit_parts.size();
		rc.cnt = tptr->parts.size();
		tptr->child_cnt[0] = 0;
		tptr->child_cnt[1] = 0;
	} else {
		auto futl = tptr->children[0].drift(stack_cnt, true, dt, abar);
		auto futr = tptr->children[1].drift(stack_cnt, false, dt, abar);
		std::array<drift_return, NCHILD> crc;
		crc[0] = futl.get();
		crc[1] = futr.get();
		std::lock_guard<mutex_type> lock(tptr->mtx);
		for (int ci = 0; ci < NCHILD; ci++) {
			auto &parts = crc[ci].parts;
			while (parts.size()) {
				const auto x = pos_to_double(parts.front().x);
				if (in_range(x, tptr->box)) {
					tptr->parts.insert(parts.front());
				} else {
					exit_parts.insert(parts.front());
				}
				parts.remove(parts.begin());
			}
		}
		rc.ndrift = exit_parts.size() + crc[0].ndrift + crc[1].ndrift;
		rc.cnt = tptr->parts.size() + crc[0].cnt + crc[1].cnt;
		rc.ekin = crc[0].ekin + crc[1].ekin;
		tptr->child_cnt[0] = crc[0].cnt;
		tptr->child_cnt[1] = crc[1].cnt;
	}
	if (tptr->level == dist_level) {
		directory->find_home(std::move(exit_parts)).get();
	}
	return rc;
}

std::atomic<int> counter(0);

int tree::place_parts(bucket &&parts) {
	std::lock_guard<mutex_type> lock(tptr->mtx);
	while (parts.size()) {
		assert(in_range(pos_to_double(parts.front().x), tptr->box));
		tptr->parts.insert(parts.front());
		parts.remove(parts.begin());
	}
	return 0;
}
