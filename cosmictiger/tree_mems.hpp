/*
 * tree_mems.hpp
 *
 *  Created on: Sep 26, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_TREE_MEMS_HPP_
#define COSMICTIGER_TREE_MEMS_HPP_


struct tree_mems {
	bucket parts;
	std::array<tree_client, NCHILD> children;
	std::array<tree_client, NSIBLING> siblings;
	tree_client parent;
	tree_client self;
	box_id_type boxid;
	std::array<std::uint64_t, NCHILD> child_cnt;
	mutex_type mtx;
	std::uint8_t leaf;

	tree_mems& operator=(const tree_mems &other) {
		parts = other.parts;
		children = other.children;
		siblings = other.siblings;
		parent = other.parent;
		self = other.self;
		boxid = other.boxid;
		child_cnt = other.child_cnt;
		leaf = other.leaf;
		return *this;
	}

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & parts;
		arc & children;
		arc & siblings;
		arc & parent;
		arc & self;
		arc & child_cnt;
		arc & boxid;
		arc & self;
	}
};



#endif /* COSMICTIGER_TREE_MEMS_HPP_ */
