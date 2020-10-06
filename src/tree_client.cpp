/*
 * tree_client.cpp
 *
 *  Created on: Sep 26, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/hpx.hpp>
#include <cosmictiger/tree.hpp>

#define MAX_STACK 9

static std::atomic<int> thread_cnt(0);
static const int max_threads = 4 * std::thread::hardware_concurrency();

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
			future = hpx::make_ready_future(hpx::async([](Function &&function, Args &&... args) {
				auto rc = function(0, std::forward<Args>(args)...);
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
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::get_parts_action>(id);
	} else {
		return hpx::make_ready_future(reinterpret_cast<tree*>(ptr)->get_parts());
	}
}

check_pair tree_client::get_child_checks() const {
	assert(hpx::get_colocation_id(id).get() == hpx::find_here());
	return reinterpret_cast<tree*>(ptr)->get_child_checks();
}

hpx::future<std::uint64_t> tree_client::prune(int stack_cnt, bool left) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::prune_action>(id, 0);
	} else {
		if (left || stack_cnt >= MAX_STACK) {
			return thread_handler<std::uint64_t>([this](int stack_cnt) {
				return reinterpret_cast<tree*>(ptr)->prune(stack_cnt);
			}, stack_cnt);
		} else {
			return hpx::make_ready_future(reinterpret_cast<tree*>(ptr)->prune(stack_cnt));
		}
	}
}
hpx::future<std::uint64_t> tree_client::grow(int stack_cnt, bool left, bucket &&parts) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::grow_action>(id, 0, std::move(parts));
	} else {
		if (left || stack_cnt >= MAX_STACK) {
			return thread_handler<std::uint64_t>([this](int stack_cnt, bucket &&parts) {
				return reinterpret_cast<tree*>(ptr)->grow(stack_cnt, std::move(parts));
			}, stack_cnt, std::move(parts));
		} else {
			return hpx::make_ready_future(reinterpret_cast<tree*>(ptr)->grow(stack_cnt, std::move(parts)));
		}
	}
}

hpx::future<multipole_return> tree_client::compute_multipoles(int stack_cnt, bool left, std::uint64_t work_id) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::compute_multipoles_action>(id, 0, work_id);
	} else {
		if (left || stack_cnt >= MAX_STACK) {
			return thread_handler<multipole_return>([this, work_id](int stack_cnt) {
				return reinterpret_cast<tree*>(ptr)->compute_multipoles(stack_cnt, work_id);
			}, stack_cnt);
		} else {
			return hpx::make_ready_future(reinterpret_cast<tree*>(ptr)->compute_multipoles(stack_cnt, work_id));
		}
	}
}

hpx::future<int> tree_client::load_balance(int stack_cnt, bool left, std::uint64_t cnt, std::uint64_t cnt2) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::load_balance_action>(id, 0, cnt, cnt2);
	} else {
		if (left || stack_cnt >= MAX_STACK) {
			return thread_handler<int>([this, cnt, cnt2](int stack_cnt) {
				return reinterpret_cast<tree*>(ptr)->load_balance(stack_cnt, cnt, cnt2);
			}, stack_cnt);
		} else {
			return hpx::make_ready_future(reinterpret_cast<tree*>(ptr)->load_balance(stack_cnt, cnt, cnt2));
		}
	}
}

hpx::future<int> tree_client::verify(int stack_cnt, bool left) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::verify_action>(id, 0);
	} else {
		if (left || stack_cnt >= MAX_STACK) {
			return thread_handler<int>([this](int stack_cnt) {
				return reinterpret_cast<tree*>(ptr)->verify(stack_cnt);
			}, stack_cnt);
		} else {
			return hpx::make_ready_future(reinterpret_cast<tree*>(ptr)->verify(stack_cnt));
		}
	}
}

hpx::future<tree_client> tree_client::migrate(hpx::id_type locality) const {

	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::migrate_action>(id, std::move(locality));
	} else {
		return hpx::make_ready_future(reinterpret_cast<tree*>(ptr)->migrate(std::move(locality)));
	}

}

hpx::future<int> tree_client::destroy(int stack_cnt) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::destroy_action>(id, 0);
	} else {
		return thread_handler<int>([this](int stack_cnt) {
			return reinterpret_cast<tree*>(ptr)->destroy(stack_cnt);
		}, stack_cnt);
	}

}
hpx::future<std::uint64_t> tree_client::drift(int stack_cnt, bool left, int step, tree_client parent, tree_client self, float dt) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::drift_action>(id, 0, step, std::move(parent), std::move(self), dt);
	} else {
		if (left || stack_cnt >= MAX_STACK) {
			return thread_handler<std::uint64_t>([this](int stack_cnt, int step, tree_client &&parent, tree_client &&self, float dt) {
				return reinterpret_cast<tree*>(ptr)->drift(stack_cnt, std::move(step), std::move(parent), std::move(self), std::move(dt));
			}, stack_cnt, std::move(step), std::move(parent), std::move(self), std::move(dt));
		} else {
			return hpx::make_ready_future(reinterpret_cast<tree*>(ptr)->drift(stack_cnt, step, std::move(parent), std::move(self), dt));
		}
	}

}

hpx::future<int> tree_client::kick_fmm(int stack_cnt, bool left, std::vector<check_item> dchecks, std::vector<check_item> echecks, expansion_src L) {
	if (hpx::get_colocation_id(id).get() != hpx::find_here()) {
		return hpx::async<tree::kick_fmm_action>(id, 0, std::move(dchecks), std::move(echecks), std::move(L));
	} else {
		if (left || stack_cnt >= MAX_STACK) {
			return thread_handler<int>(
					[this](int stack_cnt, std::vector<check_item> &&dchecks, std::vector<check_item> &&echecks, expansion_src &&L) {
						return reinterpret_cast<tree*>(ptr)->kick_fmm(stack_cnt, std::move(dchecks), std::move(echecks), std::move(L));
					}, stack_cnt, std::move(dchecks), std::move(echecks), std::move(L));
		} else {
			return hpx::make_ready_future(reinterpret_cast<tree*>(ptr)->kick_fmm(stack_cnt, std::move(dchecks), std::move(echecks), std::move(L)));
		}
	}

}

int tree_client::find_home_parent(int stack_cnt, bucket &&b) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here() || stack_cnt >= MAX_STACK) {
		auto rc = hpx::async<tree::find_home_parent_action>(id, 0, std::move(b)).get();
		return rc;
	} else {
		return reinterpret_cast<tree*>(ptr)->find_home_parent(stack_cnt, std::move(b));
	}
}

int tree_client::find_home_child(int stack_cnt, bucket &&b) const {
	if (hpx::get_colocation_id(id).get() != hpx::find_here() || stack_cnt >= MAX_STACK) {
		return hpx::async<tree::find_home_child_action>(id, 0, std::move(b)).get();
	} else {
		return reinterpret_cast<tree*>(ptr)->find_home_child(stack_cnt, std::move(b));
	}
}

std::vector<part_pos> tree_client::get_positions() const {
	assert(hpx::get_colocation_id(id).get() == hpx::find_here());
	return reinterpret_cast<tree*>(ptr)->get_positions();
}
