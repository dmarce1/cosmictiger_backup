#include <cosmictiger/hpx.hpp>
#include <cosmictiger/position.hpp>
#include <cosmictiger/tree.hpp>
#include <cosmictiger/pos_cache.hpp>

#define CACHE_WIDTH (32 * 1024)
#define CACHE_DEPTH 2

#define CACHE_NOTFOUND 0
#define CACHE_NOTREADY 1
#define CACHE_READY 2

struct cache_entry {
	tree_client id;
	bool ready;
	std::shared_ptr<std::vector<part_pos>> data_ptr;
	cache_entry() {
		id = tree_client();
		ready = false;
		data_ptr = std::make_shared<std::vector<part_pos>>();
	}
};

struct cache_line {
	std::array<std::shared_ptr<cache_entry>, CACHE_DEPTH> entries;
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
	std::vector < std::shared_ptr < cache_entry >> cache_entries(ids.size());
	std::unordered_map<int, request_type> requests;
	std::vector<bool> notfound(ids.size());
	for (int i = 0; i < ids.size(); i++) {
		const int index = gen_index(ids[i]);
		auto &line = cache[index];
		std::lock_guard<mutex_type> lock(line.mtx);
		std::shared_ptr<cache_entry> entry = nullptr;
		for (int i = 0; i < CACHE_DEPTH; i++) {
			if (line.entries[i] != nullptr) {
				if (line.entries[i]->id == ids[i]) {
					entry = line.entries[i];
					break;
				}
			}
		}
		if (entry == nullptr) {
			notfound[i] = true;
			for (int i = CACHE_DEPTH - 1; i > 0; i--) {
				line.entries[i] = line.entries[i - 1];
			}
			auto &new_entry = line.entries[0];
			new_entry = std::make_shared<cache_entry>();
			new_entry->id = ids[i];
			new_entry->ready = false;
			cache_entries[i] = new_entry;
		} else {
			notfound[i] = false;
			res[i] = entry->data_ptr;
			cache_entries[i] = entry;
		}
	}
	for (int i = 0; i < ids.size(); i++) {
		if (notfound[i]) {
			const int rank = hpx::get_locality_id_from_id(ids[i]);
			auto &entry = requests[rank];
			entry.order.push_back(i);
			entry.id.push_back(ids[i]);
		}
	}
	std::vector < hpx::future < std::vector<std::vector<part_pos>> >> futs;
	for (auto req : requests) {
		futs.push_back(hpx::async<get_remote_positions_action>(hpx_localities()[req.first], req.second.id));
	}
	int i = 0;
	for (auto req : requests) {
		auto pos = futs[i].get();
		int j = 0;
		for (auto this_pos : pos) {
			const auto resi = req.second.order[j];
			res[resi] = nullptr;
			*(cache_entries[resi]->data_ptr) = std::move(this_pos);
			cache_entries[resi]->ready = true;
			j++;
		}
		i++;
	}
	for (int i = 0; i < ids.size(); i++) {
		if (!cache_entries[i]->ready) {
			hpx::this_thread::yield();
		}
	}
	return std::move(res);
}

