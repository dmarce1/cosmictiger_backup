/*
 * tree_client.hpp
 *
 *  Created on: Sep 26, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_TREE_CLIENT_HPP_
#define COSMICTIGER_TREE_CLIENT_HPP_


class family_check;


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
	tree_client(hpx::id_type &&myid, tree* local_ptr);
	tree_client(hpx::id_type &&myid, std::uint64_t local_ptr);
	operator hpx::id_type() const {
		return id;
	}
	bool operator==( const tree_client& other ) const {
		return id == other.id;
	}
	hpx::future<int> drift(int, int,tree_client, tree_client, float dt) const;
	hpx::future<int> find_home(int, bucket&&) const;
	hpx::future<bucket> get_parts() const;
	hpx::future<int> grow(int, bucket&&) const;
	hpx::future<tree_client> migrate(hpx::id_type) const;
	hpx::future<int> load_balance(int, std::uint64_t) const;
	hpx::future<std::uint64_t> prune(int) const;
	hpx::future<int> verify(int) const;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & id;
		arc & ptr;
	}
};



#endif /* COSMICTIGER_TREE_CLIENT_HPP_ */