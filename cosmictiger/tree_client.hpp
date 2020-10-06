/*
 * tree_client.hpp
 *
 *  Created on: Sep 26, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_TREE_CLIENT_HPP_
#define COSMICTIGER_TREE_CLIENT_HPP_

#include <cosmictiger/bucket.hpp>
#include <cosmictiger/expansion.hpp>

class tree;
class check_item;
class multipole_return;
class check_info;

using check_pair = std::pair<check_item,check_item>;

class tree_client {
	hpx::id_type id;
	std::uint64_t ptr;
public:
	tree_client(tree_client&&) = default;
	tree_client(const tree_client&) = default;
	tree_client& operator=(tree_client&&) = default;
	tree_client& operator=(const tree_client&) = default;
	tree_client() {
		id = hpx::invalid_id;
		ptr = 0;
	}
	hpx::id_type get_id() const {
		return id;
	}
	std::uint64_t get_ptr() const {
		return ptr;
	}
	tree_client(hpx::id_type myid, tree *local_ptr);
	tree_client(hpx::id_type myid, std::uint64_t local_ptr);
	operator hpx::id_type() const {
		return id;
	}
	bool operator!=(const tree_client &other) const {
		return id != other.id;
	}
	bool operator==(const tree_client &other) const {
		return id == other.id;
	}
	int get_rank() const {
		return hpx::get_locality_id_from_id(id);
	}
	hpx::future<int> destroy(int) const;
	hpx::future<std::uint64_t> drift(int, bool, int, tree_client, tree_client, float dt) const;
	int find_home_parent(int, bucket&&) const;
	int find_home_child(int, bucket&&) const;
	check_pair get_child_checks() const;
	hpx::future<bucket> get_parts() const;
	hpx::future<std::uint64_t> grow(int, bool, bucket&&) const;
	hpx::future<tree_client> migrate(hpx::id_type) const;
	hpx::future<int> load_balance(int, bool left, std::uint64_t, std::uint64_t) const;
	hpx::future<multipole_return> compute_multipoles(int, bool left, std::uint64_t work_id) const;
	hpx::future<std::uint64_t> prune(int, bool) const;
	hpx::future<int> verify(int, bool) const;
	std::vector<part_pos> get_positions() const;
	check_info get_check_info() const;
	hpx::future<int> kick_fmm(int stack_cnt, bool, std::vector<check_item> dchecks, std::vector<check_item> echecks, expansion_src L);
	bool local() const {
		return hpx::get_colocation_id(id).get() == hpx::find_here();
	}
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & id;
		arc & ptr;
	}
};

#endif /* COSMICTIGER_TREE_CLIENT_HPP_ */
