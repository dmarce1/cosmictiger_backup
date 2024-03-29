#include <cosmictiger/check.hpp>
#include <cosmictiger/hpx.hpp>
#include <cosmictiger/tree.hpp>

#include <stack>
#include <unordered_map>

#define CACHE_SIZE 1024

struct cache_entry {
	bool ready;
	check_pair checks;
};

using map_type = std::unordered_map<tree_ptr, cache_entry, tree_ptr_hash>;
static map_type cache[CACHE_SIZE];
static mutex_type cache_mutex[CACHE_SIZE];
static std::stack<void*, std::vector<void*>> stack;
static mutex_type mtx;

int gen_index(const tree_ptr &id) {
	return (id.ptr >> 3) % CACHE_SIZE;
}

std::vector<check_pair> get_check_pair_remote(const std::vector<tree_ptr> &ids) {
	std::vector<check_pair> checks(ids.size());
	for (int i = 0; i < ids.size(); i++) {
		checks[i] = reinterpret_cast<tree*>(ids[i].ptr)->get_child_checks();
	}
	return std::move(checks);
}

HPX_PLAIN_ACTION (get_check_pair_remote);

hpx::future<std::vector<check_item>> get_next_checklist(const std::vector<check_item> &old) {
	std::vector<check_item> next;
	std::unordered_map<int, std::vector<tree_ptr>> get_list;
	std::vector<tree_ptr> wait_list;
	next.reserve(NCHILD * old.size());
	for (int i = 0; i < old.size(); i++) {
//		assert(old[i].opened == false);
		if (old[i].info->leaf) {
			auto j = old[i];
			j.opened = true;
			next.push_back(j);
		} else {
			const auto id = old[i].info->node;
			if (id.rank == hpx::get_locality_id()) {
				auto checks = reinterpret_cast<tree*>(id.ptr)->get_child_checks();
				next.push_back(checks.first);
				next.push_back(checks.second);
			} else {
				const int index = gen_index(id);
				std::lock_guard<mutex_type> lock(cache_mutex[index]);
				if (cache[index].find(id) == cache[index].end()) {
					get_list[id.rank].push_back(id);
					cache[index][id].ready = false;
				} else {
					const auto &entry = cache[index][id];
					if (entry.ready) {
						next.push_back(entry.checks.first);
						next.push_back(entry.checks.second);
					} else {
						wait_list.push_back(id);
					}
				}
			}
		}
	}
	if (get_list.size() || wait_list.size()) {
		return hpx::async([](decltype(next) next, decltype(get_list) get_list, decltype(wait_list) wait_list) {
			std::vector < hpx::future<std::vector<check_pair>> > futs;
			for (auto &get : get_list) {
				//			printf( "futspushback %i\n", get.first);
				futs.push_back(hpx::async<get_check_pair_remote_action>(hpx_localities()[get.first], get.second));
			}

			int i = 0;
			for (auto &get : get_list) {
				auto this_list = futs[i].get();
				//		printf( "%i ready\n", this_list.size());
				for (int j = 0; j < this_list.size(); j++) {
					next.push_back(this_list[j].first);
					next.push_back(this_list[j].second);
					const auto id = get.second[j];
					const int index = gen_index(id);
					std::lock_guard<mutex_type> lock(cache_mutex[index]);
					cache[index][id].checks = this_list[j];
					cache[index][id].ready = true;
				}
				i++;
			}

			for (auto &id : wait_list) {
				const int index = gen_index(id);
				bool done = false;
				do {
					std::unique_lock<mutex_type> lock(cache_mutex[index]);
					auto &entry = cache[index][id];
					if (entry.ready) {
						next.push_back(entry.checks.first);
						next.push_back(entry.checks.second);
						done = true;
					} else {
						lock.unlock();
						hpx::this_thread::yield();
					}
				} while (!done);
			}
			return std::move(next);
		}, std::move(next), std::move(get_list), std::move(wait_list));
	} else {
		return hpx::make_ready_future(std::move(next));
	}
}

void* check_allocate(std::size_t size) {
	auto *ptr = malloc(size);
	std::lock_guard<mutex_type> lock(mtx);
	stack.push(ptr);
	return ptr;
}

void check_cleanup();

HPX_PLAIN_ACTION (check_cleanup);

void check_cleanup() {
	std::vector<hpx::future<void>> futs;
	const auto il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const auto ir = ((hpx::get_locality_id() + 1) << 1);
	if (il < hpx_localities().size()) {
		futs.push_back(hpx::async<check_cleanup_action>(hpx_localities()[il]));
	}
	if (ir < hpx_localities().size()) {
		futs.push_back(hpx::async<check_cleanup_action>(hpx_localities()[ir]));
	}

	while (stack.size()) {
		free(stack.top());
		stack.pop();
	}
	for (int i = 0; i < CACHE_SIZE; i++) {
		cache[i].clear();
	}
	hpx::wait_all(futs.begin(), futs.end());
}

