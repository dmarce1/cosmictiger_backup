#include <cosmictiger/hpx.hpp>
#include <cosmictiger/gravity_queue.hpp>
#include <cosmictiger/pos_cache.hpp>

#include <unordered_map>
#include <unordered_set>

#define MAP_SIZE 1024

struct work_subunit1 {
	std::shared_ptr<std::vector<_4force>> f;
	std::shared_ptr<std::vector<part_pos>> x;
	std::vector<tree_client> y;
	std::vector<multi_src*> z;
	std::function<void(void)> callback;
};

struct work_subunit2 {
	std::shared_ptr<std::vector<_4force>> f;
	std::shared_ptr<std::vector<part_pos>> x;
	std::vector<std::pair<int, int>> yiter;
	std::vector<multi_src*> z;
};

struct work_unit {
	std::vector<work_subunit1> units;
	int member_count;
	int members_in;
	work_unit() {
		member_count = members_in = 0;
	}
};

struct tree_client_hash {
	std::size_t operator()(const tree_client &id) const {
		return (id.get_ptr() >> 3);
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
		futs.push_back(hpx::async < gravity_queue_retire_futures_action > (hpx_localities()[il]));
	}
	if (ir < hpx_localities().size()) {
		futs.push_back(hpx::async < gravity_queue_retire_futures_action > (hpx_localities()[ir]));
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

void gravity_queue_add_work(std::uint64_t id, std::shared_ptr<std::vector<_4force>> f, std::shared_ptr<std::vector<part_pos>> x, std::vector<tree_client> &&y,
		std::vector<multi_src*> &&z, std::function<void(void)> &&callback) {
	std::unique_lock<mutex_type> lock(mtx[id % MAP_SIZE]);
	auto *entry = &(map[id % MAP_SIZE][id / MAP_SIZE]);
	entry->members_in++;
	work_subunit1 subunit;
	subunit.f = f;
	subunit.x = x;
	subunit.y = std::move(y);
	subunit.z = std::move(z);
	subunit.callback = std::move(callback);
	entry->units.push_back(std::move(subunit));
	if (entry->members_in = entry->member_count) {
		futures.push_back(hpx::async([entry]() {
			auto unit = std::move(*entry);
			std::unordered_set<tree_client, tree_client_hash> part_requests;
			for (int i = 0; i < unit.units.size(); i++) {
				for (auto &this_y : unit.units[i].y) {
					part_requests.insert(this_y);
				}
			}
			std::vector<tree_client> part_req_vec(part_requests.begin(), part_requests.end());
			auto pos = get_positions(part_req_vec);
			std::vector < part_pos > y;
			std::vector<std::pair<int, int>> yiter;
			std::unordered_map<tree_client, std::pair<int, int>, tree_client_hash> ymap;
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
			std::vector<work_subunit2> gunits;
			for (auto &u : unit.units) {
				work_subunit2 this_subunit;
				this_subunit.f = u.f;
				this_subunit.x = u.x;
				this_subunit.z = std::move(u.z);
				std::vector<std::pair<int, int>> these_iters;
				for (int j = 0; j < u.y.size(); j++) {
					these_iters.push_back(ymap[u.y[j]]);
				}
				this_subunit.yiter = std::move(these_iters);

				gunits.push_back(std::move(this_subunit));
			}

			// CALL GRAVITY

			for (auto &u : unit.units) {
				u.callback();
			}

		}));
	}
}

