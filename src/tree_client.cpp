/*
 * tree_client.cpp
 *
 *  Created on: Sep 26, 2020
 *      Author: dmarce1
 */


#include <cosmictiger/hpx.hpp>
#include <cosmictiger/tree.hpp>



#define MAX_STACK 16


static std::atomic<int> thread_cnt(0);
static const int max_threads = 4 * std::thread::hardware_concurrency();

static bool thread_me(int stack_cnt, hpx::id_type gid) {
	bool threadme;
	if (thread_cnt++ < max_threads) {
		threadme = true;
	} else {
		thread_cnt--;
		threadme = false;
	}
	if (!threadme) {
		if (gid != hpx::find_here()) {
			thread_cnt++;
			threadme = true;
		} else if (stack_cnt >= MAX_STACK) {
			threadme = true;
		}
	}
	return threadme;
}

template<class T, class Function, class ... Args>
hpx::future<T> thread_handler(const hpx::id_type &id, const hpx::id_type loc_id, int stack_cnt, Args &&... args) {
	Function function;
	hpx::future<T> future;
	if (hpx::find_here() != loc_id) {
		future = hpx::async<Function>(id, 0, std::forward<Args>(args)...);
	} else {
		if (thread_cnt++ < max_threads) {
			future = hpx::async<Function>(id, 0, std::forward<Args>(args)...).then([](decltype(future) f) {
				thread_cnt--;
				return f.get();
			});
		} else {
			thread_cnt--;
			if (stack_cnt >= MAX_STACK) {
				future = hpx::make_ready_future(hpx::async<Function>(id, 0, std::forward<Args>(args)...).get());
			} else {
				future = hpx::make_ready_future(function(id, 0, std::forward<Args>(args)...));
			}
		}
	}
	return future;
}

tree_client::tree_client(hpx::id_type &&myid) {
	id = std::move(myid);
}

tree_client& tree_client::operator=(const hpx::id_type &myid) {
	id = std::move(myid);
	return *this;
}


hpx::future<int> tree_client::drift(int stack_cnt, int step, float dt) const {
	return thread_handler<int, tree::drift_action>(id, hpx::get_colocation_id(id).get(), stack_cnt, step, dt);
}

hpx::future<std::array<family_check, NCHILD>> tree_client::get_family_checks() const {
	return hpx::async < tree::get_family_checks_action > (id);
}

hpx::future<bucket> tree_client::get_parts() const {
	return hpx::async < tree::get_parts_action > (id);
}

hpx::future<void> tree_client::destroy() const {
	return hpx::async < tree::destroy_action > (id);
}

hpx::future<int> tree_client::find_family(int stack_cnt, tree_client parent, tree_client self, std::vector<family_check> checks) const {
	return thread_handler<int, tree::find_family_action>(id, hpx::get_colocation_id(id).get(), stack_cnt, std::move(parent), std::move(self), std::move(checks));
}

hpx::future<std::uint64_t> tree_client::grow(int stack_cnt, bucket parts) const {
	return thread_handler<std::uint64_t, tree::grow_action>(id, hpx::get_colocation_id(id).get(), stack_cnt, std::move(parts));
}

hpx::future<int> tree_client::load_balance(int stack_cnt, std::uint64_t cnt) const {
	return thread_handler<int, tree::load_balance_action>(id, hpx::get_colocation_id(id).get(), stack_cnt, cnt);
}

hpx::future<int> tree_client::verify(int stack_cnt) const {
	return thread_handler<int, tree::verify_action>(id, hpx::get_colocation_id(id).get(), stack_cnt);
}

hpx::future<int> tree_client::prune(int stack_cnt) const {
	return thread_handler<int, tree::prune_action>(id, hpx::get_colocation_id(id).get(), stack_cnt);
}

hpx::future<void> tree_client::drift_into(int sender, bucket &&parts) const {
	return hpx::async < tree::drift_into_action > (id, sender, std::move(parts));
}


hpx::future<tree_client> tree_client::migrate(hpx::id_type locality) const {
	return hpx::async<tree::migrate_action>(id, std::move(locality));
}
