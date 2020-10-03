/*
 * tree_dir_impl.hpp
 *
 *  Created on: Oct 1, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_TREE_DIR_IMPL_HPP_
#define COSMICTIGER_TREE_DIR_IMPL_HPP_


template<class A>
void tree_dir::load(A &&arc, unsigned) {
	arc & nx;
	arc & ny;
	arc & nz;
	arc & level;
	int sz;
	arc & sz;
	for (int i = 0; i < sz; i++) {
		int j;
		tree_client c;
		arc & j;
		arc & c;
		nodes[j] = c;
	}
}

template<class A>
void tree_dir::save(A &&arc, unsigned) const {
	arc & nx;
	arc & ny;
	arc & nz;
	arc & level;
	int sz = nodes.size();
	arc & sz;
	for (auto i : nodes) {
		arc & i.first;
		arc & i.second;
	}
}



#endif /* COSMICTIGER_TREE_DIR_IMPL_HPP_ */
