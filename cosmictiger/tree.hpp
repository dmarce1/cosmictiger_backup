#pragma once

#include <cosmictiger/bucket.hpp>
#include <cosmictiger/hpx.hpp>
#include <cosmictiger/range.hpp>
#include <cosmictiger/simd.hpp>

#include <array>
#include <vector>

class tree_mems;
class family_check;
class tree_client;

#include <cosmictiger/check.hpp>
#include <cosmictiger/tree_client.hpp>
#include <cosmictiger/tree_dir.hpp>

#define TREE_OVERFLOW (0x1)
#define TREE_UNDERFLOW (0x2)
#define TREE_INVALID (0x4)

std::string tree_verification_error(int rc);

struct multipole_return {
	check_info info;
	range prange;
	std::uint64_t nactive;

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & info;
		arc & nactive;
		arc & prange;
	}
};

struct fmm_params {
	float theta;
	int min_rung;
	bool stats;

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & theta;
		arc & min_rung;
		arc & stats;
	}
};

struct tree_stats {
	std::uint64_t nmig;
	std::uint64_t nnode;
	std::uint64_t nleaf;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & nmig;
		arc & nnode;
		arc & nleaf;
	}
	tree_stats() {
		nmig = nnode = nleaf = 0;
	}
	tree_stats& operator+=(const tree_stats &other) {
		nmig += other.nmig;
		nnode += other.nnode;
		nleaf += other.nleaf;
		return *this;
	}
};

void tree_set_fmm_params(fmm_params);

struct drift_in_return {
	bucket parts;
	std::uint64_t cnt;
	template<class A>
	void serialize(A&& arc, unsigned) {
		arc & parts;
		arc & cnt;
	}
};

class tree: public hpx::components::managed_component_base<tree> {
	tree_mems *tptr;
public:
	tree();
	tree(box_id_type, const range&, int);
	tree(const tree&);
	~tree();
	tree_dir build_tree_dir(tree_client) const;
	multipole_return compute_multipoles(int, std::uint64_t, std::uint64_t);
	void create_children();
	int destroy(int);
	std::uint64_t drift(int, int, tree_client, tree_client, float dt);
	drift_in_return drift_in(int, float dt);
	int drift_out(int, bucket&&, std::uint64_t);
	int find_home_parent(int, bucket);
	int find_home_child(int, bucket);
	check_pair get_child_checks() const;
	check_info get_check_info() const;
	std::vector<part_pos> get_positions() const;
	bucket get_parts();
	std::uint64_t get_ptr();
	std::uint64_t grow(int, bucket&&, bool);
	int kick_fmm(int, std::vector<check_item>&&, std::vector<check_item>&&, expansion_src&&);
	tree_stats load_balance(int stack_cnt, std::uint64_t index, std::uint64_t total);
	std::uint64_t local_load_balance(std::uint64_t, std::uint64_t);
	tree_client migrate(hpx::id_type);
	std::uint64_t prune(int);
	int verify(int) const;
	std::size_t size() const;
	std::vector<bool> checks_far(const std::vector<check_item>&, bool ewald);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,build_tree_dir);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,destroy);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,drift_in);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,drift_out);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,drift);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,find_home_parent);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,find_home_child);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,kick_fmm);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,grow);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,load_balance);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,compute_multipoles);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,prune);
	/**/HPX_DEFINE_COMPONENT_ACTION(tree,verify);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_child_checks);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_check_info);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_ptr);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_parts);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,get_positions);
	/**/HPX_DEFINE_COMPONENT_DIRECT_ACTION(tree,migrate);
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & *tptr;
	}

};

void tree_insert_parts(bucket &&parts);

#include <cosmictiger/tree_dir_impl.hpp>
#include <cosmictiger/tree_mems.hpp>
