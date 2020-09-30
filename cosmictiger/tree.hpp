#pragma once

#include <cosmictiger/bucket.hpp>
#include <cosmictiger/hpx.hpp>
#include <cosmictiger/range.hpp>

#include <array>
#include <vector>

class tree_mems;
class family_check;
class tree_client;
class tree_dir;

#define TREE_OVERFLOW (0x1)
#define TREE_UNDERFLOW (0x2)
#define TREE_INVALID (0x4)

std::string tree_verification_error(int rc);

class tree: public hpx::components::managed_component_base<tree> {
	tree_mems *tptr;

public:
	tree();
	tree(const range&, int);
	tree(const tree&);
	~tree();
//	tree_dir build_tree_dir(tree_client) const;
	void create_children();
	int destroy(int);
	std::uint64_t drift(int, int, tree_client, tree_client, float dt);
	int find_home_parent(int, bucket);
	int find_home_child(int, bucket);
	bucket get_parts();
	std::uint64_t get_ptr();
	std::uint64_t grow(int, bucket&&);
	int load_balance(int, std::uint64_t, std::uint64_t);
	tree_client migrate(hpx::id_type);
	std::uint64_t prune(int);
	int verify(int) const;
	std::size_t size() const;
//	/**/HPX_DEFINE_COMPONENT_ACTION(tree,build_tree_dir);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,destroy);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,drift);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,find_home_parent);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,find_home_child);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,grow);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,load_balance);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,prune);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,verify);
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
