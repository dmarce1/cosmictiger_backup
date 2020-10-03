#include <cosmictiger/hpx.hpp>
#include <cosmictiger/future_data.hpp>
#include <cosmictiger/position.hpp>
#include <cosmictiger/tree.hpp>

#define CACHE_WIDTH (32 * 1024)
#define CACHE_DEPTH 2

struct cache_entry {
	tree_client id;
	shared_future_data<std::vector<part_pos>> fut;
};

struct cache_line {
	std::array<cache_entry, CACHE_DEPTH> entries;
	mutex_type mtx;
};

static std::array<cache_line, CACHE_WIDTH> cache;

static std::size_t gen_index(tree_client id) {
	return (id.get_ptr() >> 3) % CACHE_WIDTH;
}

future_data<const std::vector<part_pos>*> pos_cache_read(tree_client id) {
	future_data<const std::vector<part_pos>*> fut;
	const auto index = gen_index(id);
	auto& line = cache[index];
	bool found = false;
	std::unique_lock<mutex_type> lock(line.mtx);
	for (int i = 0; i < CACHE_DEPTH; i++) {
		auto& entry = line.entries[i];
		if (entry.id == id) {
			found = true;
			fut.set(&entry.fut.get());
			break;
		}
	}
	if (!found) {
		for (int i = CACHE_DEPTH - 1; i > 0; i--) {
			line.entries[i] = line.entries[i - 1];
		}
		auto& entry = line.entries[0];
		if (id.local()) {
			entry.fut.set(reinterpret_cast<tree*>(id.get_ptr())->get_positions());
			fut.set(&entry.fut.get());
		} else {
			hpx::lcos::local::promise < hpx::future<std::vector<part_pos>> > promise;
			auto this_fut = promise.get_future();
			entry.fut.set(this_fut.then([](decltype(this_fut) f) {
				return f.get().get();
			}));
			fut.set(&entry.fut.get());
			lock.unlock();
			promise.set_value(hpx::async < tree::get_positions_action > (id.get_id()));
		}
	}
	return fut;
}


