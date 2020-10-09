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

struct tree_ptr {
	int rank;
	std::uint64_t ptr;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & rank;
		arc & ptr;
	}
	inline bool operator==(tree_ptr other) const {
		return rank == other.rank && ptr == other.ptr;
	}
};

struct tree_ptr_hash {
	std::size_t operator()(tree_ptr ptr) const {
		return ptr.ptr ^ ptr.rank;
	}
};

struct check_info {
	multi_src *multi;
	tree_ptr node;
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
		arc & opened;
		arc & info->leaf;
		arc & info->node;
	}
	template<class A>
	void save(A &&arc, unsigned) const {
		arc & *(info->multi);
		arc & opened;
		arc & info->leaf;
		arc & info->node;
	}

};

using check_pair = std::pair<check_item, check_item>;

hpx::future<std::vector<check_item>> get_next_checklist(const std::vector<check_item> &old);

#endif /* COSMICTIGER_CHECK_ITEM_HPP_ */
