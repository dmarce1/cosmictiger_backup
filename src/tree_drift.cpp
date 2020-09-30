/*
 * tree_drift.cpp
 *
 *  Created on: Sep 30, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/tree.hpp>

std::uint64_t tree::drift(int stack_cnt, int step, tree_client parent, tree_client self, float dt) {
	tptr->parent = parent;
	std::uint64_t drifted = 0;
	if (tptr->leaf) {
		while( tptr->lock++ != 0 ) {
			tptr->lock--;
		}
		bucket exit_parts;
		auto i = tptr->parts.begin();
		while (i != tptr->parts.end()) {
			if (step % 2 == i->step) {
				auto x = pos_to_double(i->x);
				const auto v = i->v;
				x += v * dt;
				i->x = double_to_pos(x);
				i->step++;
				if (!in_range(pos_to_double(i->x), tptr->box)) {
					exit_parts.insert(*i);
					drifted++;
					i = tptr->parts.remove(i);
				} else {
					i++;
				}
			} else {
				i++;
			}
		}
		tptr->lock--;
		if (exit_parts.size()) {
			parent.find_home_parent(stack_cnt, std::move(exit_parts));
		}
	} else {
		auto futl = tptr->children[0].drift(stack_cnt, true, step, self, tptr->children[0], dt);
		auto futr = tptr->children[1].drift(stack_cnt, false, step, self, tptr->children[1], dt);
		drifted += futl.get();
		drifted += futr.get();
	}
	return drifted;
}

int tree::find_home_child(int stack_cnt, bucket parts) {
	if (tptr->leaf) {
		while( tptr->lock++ != 0 ) {
			tptr->lock--;
		}
		while (parts.size()) {
			assert(in_range(pos_to_double(parts.front().x), tptr->box));
			tptr->parts.insert(parts.front());
			parts.remove(parts.begin());
		}
		tptr->lock--;
	} else {
		bucket l_parts;
		bucket r_parts;
		const auto dim = tptr->level % NDIM;
		const auto midx = 0.5 * (tptr->box.min[dim] + tptr->box.max[dim]);
		while (parts.size()) {
			auto &p = parts.front();
			const auto x = pos_to_double(p.x);
			if (x[dim] >= midx) {
				r_parts.insert(p);
			} else {
				l_parts.insert(p);
			}
			parts.remove(parts.begin());
		}
		if (l_parts.size()) {
			assert(tptr->children[0] != tree_client());
			tptr->children[0].find_home_child(stack_cnt, std::move(l_parts));
		}
		if (r_parts.size()) {
			assert(tptr->children[1] != tree_client());
			tptr->children[1].find_home_child(stack_cnt, std::move(r_parts));
		}
	}
	return 0;
}

int tree::find_home_parent(int stack_cnt, bucket parts) {
	bucket p_parts;
	bucket l_parts;
	bucket r_parts;
	const auto dim = tptr->level % NDIM;
	const auto midx = 0.5 * (tptr->box.min[dim] + tptr->box.max[dim]);
	while (parts.size()) {
		auto &p = parts.front();
		const auto x = pos_to_double(p.x);
		if (!in_range(x, tptr->box)) {
			p_parts.insert(p);
		} else if (x[dim] >= midx) {
			r_parts.insert(p);
		} else {
			l_parts.insert(p);
		}
		parts.remove(parts.begin());
	}
	if (p_parts.size()) {
		assert(tptr->parent != tree_client());
		tptr->parent.find_home_parent(stack_cnt, std::move(p_parts));
	}
	if (l_parts.size()) {
		assert(tptr->children[0] != tree_client());
		tptr->children[0].find_home_child(stack_cnt, std::move(l_parts));
	}
	if (r_parts.size()) {
		assert(tptr->children[1] != tree_client());
		tptr->children[1].find_home_child(stack_cnt, std::move(r_parts));
	}
	return 0;
}
