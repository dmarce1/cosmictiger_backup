#include <cosmictiger/hpx.hpp>
#include <cosmictiger/gravity_cuda.hpp>
#include <cosmictiger/gravity_queue.hpp>
#include <cosmictiger/pos_cache.hpp>
#include <cosmictiger/pinned_vector.hpp>

#include <unordered_map>
#include <unordered_set>

#define MAP_SIZE 1024

struct work_subunit {
	std::shared_ptr<std::vector<_4force>> f;
	std::shared_ptr<std::vector<part_pos>> x;
	std::vector<tree_ptr> y;
	std::vector<const multi_src*> z;
	std::function<void(void)> callback;
};

struct work_unit {
	std::vector<work_subunit> units;
	int member_count;
	int members_in;
	work_unit() {
		member_count = members_in = 0;
	}
};

static std::atomic<std::uint64_t> next_index(0);
static std::unordered_map<std::uint64_t, work_unit> map[MAP_SIZE];
static mutex_type mtx[MAP_SIZE];
static mutex_type fut_mtx;
static std::vector<hpx::future<void>> futures;

std::uint64_t gravity_queue_genid() {
	return next_index++ * std::uint64_t(hpx_localities().size()) + std::uint64_t(hpx::get_locality_id());
}

void gravity_queue_checkin(std::uint64_t id) {
	std::lock_guard<mutex_type> lock(mtx[id % MAP_SIZE]);
	map[id % MAP_SIZE][id / MAP_SIZE].member_count++;
}

HPX_PLAIN_ACTION (gravity_queue_retire_futures);

void gravity_queue_retire_futures() {
	std::vector<hpx::future<void>> futs;
	const auto il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const auto ir = ((hpx::get_locality_id() + 1) << 1);
	if (il < hpx_localities().size()) {
		futs.push_back(hpx::async<gravity_queue_retire_futures_action>(hpx_localities()[il]));
	}
	if (ir < hpx_localities().size()) {
		futs.push_back(hpx::async<gravity_queue_retire_futures_action>(hpx_localities()[ir]));
	}
	for (auto &f : futures) {
		f.get();
	}
	for (int i = 0; i < MAP_SIZE; i++) {
		map[i].clear();
	}
	futures.clear();
	hpx::wait_all(futs.begin(), futs.end());
}

void gravity_queue_add_work(std::uint64_t id, std::shared_ptr<std::vector<_4force>> f, std::shared_ptr<std::vector<part_pos>> x, std::vector<tree_ptr> &&y,
		std::vector<const multi_src*> &&z, std::function<void(void)> &&callback, bool do_phi) {
	std::unique_lock<mutex_type> lock(mtx[id % MAP_SIZE]);
	auto *entry = &(map[id % MAP_SIZE][id / MAP_SIZE]);
	entry->members_in++;
	work_subunit subunit;
	subunit.f = f;
	subunit.x = x;
	subunit.y = std::move(y);
	subunit.z = std::move(z);
	subunit.callback = std::move(callback);
	entry->units.push_back(std::move(subunit));
	assert(entry->members_in <= entry->member_count);
	if (entry->members_in == entry->member_count) {
//		auto func = hpx::async([entry, do_phi]() {
		auto unit = std::move(*entry);
		std::unordered_set<tree_ptr, tree_ptr_hash> part_requests;
		for (int i = 0; i < unit.units.size(); i++) {
			for (auto &this_y : unit.units[i].y) {
				part_requests.insert(this_y);
			}
		}
		std::vector<tree_ptr> part_req_vec(part_requests.begin(), part_requests.end());
		auto pos = get_positions(part_req_vec);
		pinned_vector<part_pos> y;
		y.reserve(opts.bucket_size * pos.size());
		std::vector<std::pair<int, int>> yiter;
		static thread_local std::unordered_map<tree_ptr, std::pair<int, int>, tree_ptr_hash> ymap;
		ymap.clear();
		int index = 0;
		int i = 0;
		for (const auto &this_pos : pos) {
			for (const auto &p : *this_pos) {
				y.push_back(p);
			}
			ymap[part_req_vec[i]] = std::make_pair(index, index + this_pos->size());
			index += this_pos->size();
			i++;
		}
		std::vector<cuda_work_unit> gunits;
		for (auto &u : unit.units) {
			cuda_work_unit this_subunit;
			this_subunit.fptr = u.f;
			this_subunit.xptr = u.x;
			this_subunit.z = std::move(u.z);
			std::vector<std::pair<int, int>> these_iters;
			these_iters.reserve(u.y.size());
			for (int j = 0; j < u.y.size(); j++) {
				these_iters.push_back(ymap[u.y[j]]);
			}
			for (auto &iter : these_iters) {
				while (iter.first != iter.second) {
					std::pair<int, int> this_iter;
					this_iter.first = iter.first;
					this_iter.second = std::min(iter.second, this_iter.first + SYNCRATE);
					iter.first = this_iter.second;
					this_subunit.yiters.push_back(this_iter);
				}

			}
			gunits.push_back(std::move(this_subunit));
		}

		gravity_PP_direct_cuda(std::move(gunits), y, do_phi);

		for (auto &u : unit.units) {
			u.callback();
		}

//		});
//		func.get();
//		std::lock_guard < mutex_type > lock(fut_mtx);
//		futures.push_back(std::move(func));
	}
}

