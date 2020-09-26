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
public:
	/**/DEFAULT_CLASS_MEMBERS(tree_client)
	;
	tree_client(hpx::id_type &&myid);
	tree_client& operator=(const hpx::id_type &myid);
	operator hpx::id_type() const {
		return id;
	}
	hpx::future<void> destroy() const;
	hpx::future<int> drift(int, int, float dt) const;
	hpx::future<void> drift_into(int, bucket&&) const;
	hpx::future<int> find_family(int, tree_client, tree_client, std::vector<family_check>) const;
	hpx::future<std::array<family_check, NCHILD>> get_family_checks() const;
	hpx::future<bucket> get_parts() const;
	hpx::future<std::uint64_t> grow(int, bucket&&) const;
	hpx::future<tree_client> migrate(hpx::id_type) const;
	hpx::future<int> load_balance(int, std::uint64_t) const;
	hpx::future<int> prune(int) const;
	hpx::future<int> verify(int) const;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & id;
	}
};



struct family_check {
	box_id_type boxid;
	tree_client node;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & boxid;
		arc & node;
	}
};


#endif /* COSMICTIGER_TREE_CLIENT_HPP_ */
