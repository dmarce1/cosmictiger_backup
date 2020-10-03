/*
 * check_item.hpp
 *
 *  Created on: Oct 2, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_CHECK_ITEM_HPP_
#define COSMICTIGER_CHECK_ITEM_HPP_

#include <cosmictiger/multipole.hpp>
#include <cosmictiger/tree_client.hpp>

void* check_allocate(std::size_t);
void check_cleanup();

struct check_info {
	multi_src *multi;
	tree_client node;
	bool leaf;

	HPX_SERIALIZATION_SPLIT_MEMBER();
	template<class A>
	void load(A &&arc, unsigned) {
		bool null;
		arc & null;
		if (!null) {
			multi = (multi_src*) check_allocate(sizeof(multi_src));
			arc & *multi;
		} else {
			multi = nullptr;
		}
		arc & node;
		arc & leaf;
	}

	template<class A>
	void save(A &&arc, unsigned) const {
		bool null = multi == nullptr;
		arc & null;
		if (!null) {
			arc & *multi;
		}
		arc & node;
		arc & leaf;
	}

};

struct check_item {
	check_info *info;
	bool opened;

	HPX_SERIALIZATION_SPLIT_MEMBER();
	template<class A>
	void load(A &&arc, unsigned) {
		auto *ptr = check_allocate(sizeof(check_info) + sizeof(multi_src));
		info = (check_info*) ptr;
		info->multi = (multi_src*) ((void*) ptr + sizeof(check_info));
		arc & *(info->multi);
		arc & info->node;
		arc & info->leaf;
		arc & opened;
	}
	template<class A>
	void save(A &&arc, unsigned) const {
		arc & *(info->multi);
		arc & info->node;
		arc & info->leaf;
		arc & opened;
	}

};

using check_pair = std::pair<check_item, check_item>;

#endif /* COSMICTIGER_CHECK_ITEM_HPP_ */
