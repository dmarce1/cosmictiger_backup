/*
 * tree_client.cpp
 *
 *  Created on: Sep 26, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/hpx.hpp>
#include <cosmictiger/tree.hpp>

#define MAX_STACK 18

static std::atomic<int> thread_cnt(0);
static const int max_threads = 4 * std::thread::hardware_concurrency();
static std::atomic<int> stack_threads(0);

template<class T, class Function, class ... Args>
hpx::future<T> thread_handler(Function &&function, int stack_cnt, Args &&... args) {
	hpx::future<T> future;
	stack_cnt++;
	if (thread_cnt++ < max_threads) {
		future = hpx::async([](Function &&function, Args &&... args) {
			auto rc = function(0, std::forward<Args>(args)...);
			thread_cnt--;
			return rc;
		}, std::move(function), std::forward<Args>(args)...);
	} else {
		thread_cnt--;
		if (stack_cnt >= MAX_STACK) {
			stack_threads++;
			future = hpx::make_ready_future(hpx::async([](Function &&function, Args &&... args) {
				 auto rc = function(0, std::forward<Args>(args)...);
				 stack_threads--;
	//			 printf( "%i\n", (int) stack_threads);
				 return rc;
			}, std::move(function), std::forward<Args>(args)...).get());
		} else {
			future = hpx::make_ready_future(function(stack_cnt, std::forward<Args>(args)...));
		}
	}
	return future;
}

tree_client::tree_client(hpx::id_type myid, tree *local_ptr) {
	ptr = reinterpret_cast<std::uint64_t>(local_ptr);
	id = std::move(myid);
}

tree_client::tree_client(hpx::id_type myid, std::uint64_t local_ptr) {
	ptr = local_ptr;
	id = std::move(myid);
}

hpx::future<bucket> tree_client::get_parts() const {
	return hpx::async < tree::get_parts_action > (id);
}

hpx::future<std::uint64_t> tree_client::prune(int stack_cnt) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async < tree::prune_action > (id, 0);
	} else {
		return thread_handler<std::uint64_t>([this](int stack_cnt) {
			return reinterpret_cast<tree*>(ptr)->prune(stack_cnt);
		}, stack_cnt);
	}
}
hpx::future<std::uint64_t> tree_client::grow(int stack_cnt, bucket &&parts) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async < tree::grow_action > (id, 0, std::move(parts));
	} else {
		return thread_handler<std::uint64_t>([this](int stack_cnt, bucket &&parts) {
			return reinterpret_cast<tree*>(ptr)->grow(stack_cnt, std::move(parts));
		}, stack_cnt, std::move(parts));
	}
}

hpx::future<int> tree_client::load_balance(int stack_cnt, std::uint64_t cnt) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async < tree::load_balance_action > (id, 0, cnt);
	} else {
		return thread_handler<int>([this, cnt](int stack_cnt) {
			return reinterpret_cast<tree*>(ptr)->load_balance(stack_cnt, cnt);
		}, stack_cnt);
	}
}

hpx::future<int> tree_client::verify(int stack_cnt) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async < tree::verify_action > (id, 0);
	} else {
		return thread_handler<int>([this](int stack_cnt) {
			return reinterpret_cast<tree*>(ptr)->verify(stack_cnt);
		}, stack_cnt);
	}
}

hpx::future<tree_client> tree_client::migrate(hpx::id_type locality) const {
	return hpx::async < tree::migrate_action > (id, std::move(locality));
}

hpx::future<int> tree_client::destroy(int stack_cnt) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async < tree::destroy_action > (id, 0);
	} else {
		return thread_handler<int>([this](int stack_cnt) {
			return reinterpret_cast<tree*>(ptr)->destroy(stack_cnt);
		}, stack_cnt);
	}

}
hpx::future<std::uint64_t> tree_client::drift(int stack_cnt, int step, tree_client parent, tree_client self, float dt) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async < tree::drift_action > (id, 0, step, std::move(parent), std::move(self), dt);
	} else {
		return thread_handler<std::uint64_t>([this, parent, self, dt, step](int stack_cnt) {
			return reinterpret_cast<tree*>(ptr)->drift(stack_cnt, step, std::move(parent), std::move(self), dt);
		}, stack_cnt);
	}

}

hpx::future<int> tree_client::find_home(int stack_cnt, bucket &&b) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async < tree::find_home_action > (id, 0, std::move(b));
	} else {
		return thread_handler<int>([this](int stack_cnt, bucket &&b) {
			return reinterpret_cast<tree*>(ptr)->find_home(stack_cnt, std::move(b));
		}, stack_cnt, std::move(b));
	}

}

