#pragma once

#include <cosmictiger/bucket.hpp>
#include <cosmictiger/hpx.hpp>
#include <cosmictiger/range.hpp>

#include <array>
#include <vector>

class tree_mems;
class family_check;
class tree_client;


#define TREE_OVERFLOW (0x1)
#define TREE_UNDERFLOW (0x2)
#define TREE_INVALID (0x4)

std::string tree_verification_error(int rc);

class tree: public hpx::components::managed_component_base<tree> {
	tree_mems *tptr;

public:
	tree();
	tree(box_id_type id);
	tree(const tree&);
	~tree();
	void create_children();
	int drift(int, int, float dt);
	void drift_into(int, bucket);
	int find_family(int, tree_client, tree_client, std::vector<family_check>);
	void destroy();
	bucket get_parts();
	std::uint64_t get_ptr();
	std::array<family_check, NCHILD> get_family_checks() const;
	std::uint64_t grow(int, bucket);
	int load_balance(int, std::uint64_t);
	tree_client migrate(hpx::id_type);
	int prune(int);
	int verify(int) const;
	std::size_t size() const;
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,destroy);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,drift);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,drift_into);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,grow);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,load_balance);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,prune);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,verify);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,find_family);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_family_checks);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_ptr);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_parts);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,migrate);
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & *tptr;
	}

};

#include <cosmictiger/tree_client.hpp>
#include <cosmictiger/tree_mems.hpp>
