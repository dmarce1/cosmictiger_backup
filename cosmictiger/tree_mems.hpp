/*
 * tree_mems.hpp
 *
 *  Created on: Sep 26, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_TREE_MEMS_HPP_
#define COSMICTIGER_TREE_MEMS_HPP_

#include <cosmictiger/multipole.hpp>
#include <cosmictiger/range.hpp>

struct tree_mems {
	multi_src multi;
	bucket parts;
	range box;
	std::array<tree_client, NCHILD> children;
	std::array<check_info, NCHILD> child_info;
	tree_client parent;
	std::array<std::uint64_t, NCHILD> child_cnt;
	box_id_type id;
	std::uint64_t nactive;
	std::uint64_t work_id;
	mutex_type mtx;
	std::uint8_t leaf;
	std::uint8_t level;

	tree_mems() {
	}

	tree_mems& operator=(const tree_mems &other) {
		multi = other.multi;
		parts.clear();
		for (auto i = other.parts.cbegin(); i != other.parts.cend(); i++) {
			parts.insert(*i);
		}
		children = other.children;
		box = other.box;
		child_info = other.child_info;
		parent = other.parent;
		child_cnt = other.child_cnt;
		id = other.id;
		nactive = other.nactive;
		work_id = other.work_id;
		nactive = other.nactive;
		leaf = other.leaf;
		level = other.level;
		return *this;
	}

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & children;
		arc & work_id;
		arc & nactive;
		arc & multi;
		arc & parts;
		arc & box;
		arc & child_info;
		arc & parent;
		arc & child_cnt;
		arc & id;
		arc & leaf;
		arc & level;
	}
};

#endif /* COSMICTIGER_TREE_MEMS_HPP_ */
