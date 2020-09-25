#pragma once

#include <cosmictiger/bucket.hpp>
#include <cosmictiger/hpx.hpp>
#include <cosmictiger/range.hpp>

#include <array>
#include <vector>

class tree_mems;


class tree: public hpx::components::managed_component_base<tree> {
	tree_mems *tptr;

public:
	tree();
	tree(box_id_type id);
	~tree();
	void create_children();
	void add_parts(std::vector<particle>&&, int);
	HPX_DEFINE_COMPONENT_ACTION(tree,add_parts);
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & *tptr;
	}

};

class tree_client {
	hpx::id_type id;
public:
	/**/DEFAULT_CLASS_MEMBERS(tree_client);
	tree_client(hpx::id_type &&myid);
	hpx::future<void> add_parts(std::vector<particle>&&, int);
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & id;
	}
};


struct tree_mems {
	bucket parts;
	std::array<tree_client, NCHILD> children;
	std::array<tree_client, NSIBLING> siblings;
	box_id_type boxid;
	std::uint8_t leaf;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & parts;
		arc & children;
		arc & siblings;
		arc & boxid;
	}
};
