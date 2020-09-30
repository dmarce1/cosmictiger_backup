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
	range box;
	std::array<tree_client, NCHILD> children;
	tree_client parent;
	std::array<std::uint64_t, NCHILD> child_cnt;
	std::atomic<int> lock;
	std::uint8_t leaf;
	std::uint8_t level;

	tree_mems() : lock(0) {
	}

	tree_mems& operator=(const tree_mems& other) {
		parts.clear();
		for( auto i = other.parts.cbegin(); i != other.parts.cend(); i++) {
			parts.insert(*i);
		}
		children = other.children;
		parent = other.parent;
		box = other.box;
		level = other.level;
		child_cnt = other.child_cnt;
		leaf = other.leaf;
		return *this;
	}

	template<class A>
	void serialize(A &&arc, unsigned) {
	//	printf( "Serialize\n");
		arc & parts;
		arc & box;
		arc & children;
		arc & parent;
		arc & child_cnt;
		arc & box;
		arc & leaf;
		arc & level;
	}
};



#endif /* COSMICTIGER_TREE_MEMS_HPP_ */
