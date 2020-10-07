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
#include <cosmictiger/tree_client.hpp>

#include <unordered_map>

class tree_dir {
	std::unordered_map<int, tree_client> nodes;
	int nx;
	int ny;
	int nz;
	int level;
	int index(int x, int y, int z) const;
	int index(const vect<double> &x) const;
	mutex_type mtx;
public:
	tree_dir();

	tree_dir(tree_dir &&other) {
		*this = std::move(other);
	}

	tree_dir& operator=(tree_dir &&other) {
		nodes = std::move(other.nodes);
		nx = other.nx;
		ny = other.ny;
		nz = other.nz;
		return *this;
	}
	tree_dir(const tree_dir &other) {
		*this = other;
	}

	tree_dir& operator=(const tree_dir &other) {
		nodes = other.nodes;
		nx = other.nx;
		ny = other.ny;
		nz = other.nz;
		return *this;
	}

	int size() const {
		return nodes.size();
	}

	void add_tree_client(const tree_client&, const range &box);

	tree_dir& merge(const tree_dir &other);

	hpx::future<void> find_home(bucket &&p);
	HPX_SERIALIZATION_SPLIT_MEMBER();
	template<class A>
	void load(A &&arc, unsigned);

	template<class A>
	void save(A &&arc, unsigned) const;
};

#endif /* COSMICTIGER_TREE_DIR_HPP_ */
