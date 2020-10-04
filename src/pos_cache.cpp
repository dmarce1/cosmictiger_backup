#include <cosmictiger/hpx.hpp>
#include <cosmictiger/future_data.hpp>
#include <cosmictiger/position.hpp>
#include <cosmictiger/tree.hpp>

#define CACHE_WIDTH (32 * 1024)
#define CACHE_DEPTH 2

#define CACHE_NOTFOUND 0
#define CACHE_NOTREADY 1
#define CACHE_READY 2

struct cache_entry {
	tree_client id;
	bool ready;
	std::shared_ptr<std::vector<part_pos>> data_ptr;
};

struct cache_line {
	std::array<cache_entry, CACHE_DEPTH> entries;
	mutex_type mtx;
};

static std::array<cache_line, CACHE_WIDTH> cache;

static std::size_t gen_index(tree_client id) {
	return (id.get_ptr() >> 3) % CACHE_WIDTH;
}

struct request_type {
	std::vector<int> order;
	std::vector<tree_client> id;
};

std::vector<std::vector<part_pos>> get_remote_positions(const std::vector<tree_client> ids);

HPX_PLAIN_ACTION (get_remote_positions);

std::vector<std::vector<part_pos>> get_remote_positions(const std::vector<tree_client> ids) {
	std::vector < std::vector < part_pos >> pos(ids.size());
	for (int i = 0; i < ids.size(); i++) {
		pos.push_back(ids[i].get_positions());
	}
	return std::move(pos);

}

std::vector<std::shared_ptr<std::vector<part_pos>>> get_positions(const std::vector<tree_client> &ids) {
	std::vector < std::shared_ptr<std::vector<part_pos>> > res(ids.size());
	std::vector<int> status(ids.size());
	std::unordered_map<int, request_type> requests;
	for (int i = 0; i < ids.size(); i++) {
		const int index = gen_index(ids[i]);
		auto &line = cache[index];
		std::lock_guard<mutex_type> lock(line.mtx);
		cache_entry *entry = nullptr;
		for (int i = 0; i < CACHE_DEPTH; i++) {
			if (line.entries[i].id == ids[i]) {
				entry = &(line.entries[i]);
				break;
			}
		}
		if (entry == nullptr) {
			status[i] = CACHE_NOTFOUND;
			for (int i = CACHE_DEPTH - 1; i > 0; i--) {
				line.entries[i] = line.entries[i - 1];
			}
			auto &entry = line.entries[0];
			entry.id = ids[i];
			entry.ready = false;
		} else if (entry->ready) {
			status[i] = CACHE_READY;
			res[i] = entry->data_ptr;
		} else {
			status[i] = CACHE_NOTREADY;
		}
	}
	for (int i = 0; i < ids.size(); i++) {
		if (status[i] == CACHE_NOTFOUND) {
			const int rank = hpx::get_locality_id_from_id(ids[i]);
			auto &entry = requests[rank];
			entry.order.push_back(i);
			entry.id.push_back(ids[i]);
		}
	}
	std::vector < hpx::future < std::vector<std::vector<part_pos>> >> futs;
	for (auto req : requests) {
		futs.push_back(hpx::async < get_remote_positions_action > (hpx_localities()[req.first], req.second.id));
	}
	int i = 0;
	for (auto req : requests) {
		auto pos = futs[i].get();
		int j = 0;
		for (auto this_pos : pos) {
			auto ptr = std::make_shared < std::vector < part_pos >> (std::move(this_pos));
			const auto resi = req.second.order[j];
			res[resi] = ptr;
			const int index = gen_index(ids[i]);
			auto &line = cache[index];
			std::lock_guard<mutex_type> lock(line.mtx);
			for (int i = 0; i < CACHE_DEPTH; i++) {
				if (line.entries[i].id == ids[resi]) {
					line.entries[i].data_ptr = ptr;
					line.entries[i].ready = true;
				}
			}
			j++;
		}
		i++;
	}
	return std::move(res);
}

