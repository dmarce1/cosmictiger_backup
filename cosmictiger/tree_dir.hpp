/*
 * tree_dir.hpp
 *
 *  Created on: Sep 29, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_TREE_DIR_HPP_
#define COSMICTIGER_TREE_DIR_HPP_

#include <cosmictiger/defs.hpp>
#include <cosmictiger/range.hpp>
#include <cosmictiger/tree.hpp>

#include <unordered_map>

class tree_dir {
	std::unordered_map<int, tree_client> nodes;
	int nx;
	int ny;
	int nz;
	int level;
	int index(int x, int y, int z) const;
	int index(const vect<double> &x) const;
public:
	tree_dir();

	void add_tree_client(const tree_client&, box_id_type id);

	tree_dir& merge(const tree_dir &other);

	hpx::future<void> find_home(int, bucket &&p);

	HPX_SERIALIZATION_SPLIT_MEMBER();
	template<class A>
	void load(A &&arc, unsigned) {
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
	void save(A &&arc, unsigned) const {
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

};

#endif /* COSMICTIGER_TREE_DIR_HPP_ */
