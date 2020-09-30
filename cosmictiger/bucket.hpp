/*
 * bucket.hpp
 *
 *  Created on: Sep 25, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_BUCKET_HPP_
#define COSMICTIGER_BUCKET_HPP_

#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/particle.hpp>

#include <cassert>

class bucket;

#define CUP_SIZE 16

struct cup {
	std::array<particle, CUP_SIZE> data;
	cup *next;
	cup *prev;
};

cup* bucket_cup_alloc();
void bucket_cup_free(cup *ptr);

struct bucket {
	class iterator {
	public:
		cup *ptr;
		int off;
	public:
		/**/DEFAULT_CLASS_MEMBERS(iterator)
		;
		inline particle operator*() const {
			assert(ptr);
			return ptr->data[off];
		}
		inline particle& operator*() {
			assert(ptr);
			return ptr->data[off];
		}
		inline particle* operator->() {
			return &(ptr->data[off]);
		}
		inline const particle* operator->() const {
			return &(ptr->data[off]);
		}
		inline bool operator==(const iterator &other) const {
			return (ptr == other.ptr) && (off == other.off);
		}
		inline bool operator!=(const iterator &other) const {
			return !(*this == other);
		}
		inline iterator operator++() {
			assert(ptr);
			off++;
			if (off == CUP_SIZE && ptr->next) {
				off = 0;
				ptr = ptr->next;
			}
			return *this;
		}
		inline iterator operator++(int) {
			return ++(*this);
		}
		friend class bucket;
	};
	class const_iterator {
		iterator iter;
	public:
		/**/DEFAULT_CLASS_MEMBERS(const_iterator)
		;
		inline const_iterator(iterator i) {
			iter = i;
		}
		inline particle operator*() const {
			return *iter;
		}
		inline const particle* operator->() const {
			return iter.operator->();
		}
		inline bool operator==(const const_iterator &other) const {
			return (iter == other.iter);
		}
		inline bool operator!=(const const_iterator &other) const {
			return !(*this == other);
		}
		inline const_iterator operator++() {
			iter++;
			return *this;
		}
		inline const_iterator operator++(int) {
			return ++(*this);
		}
		friend class bucket;
	};
private:
	iterator start;
	iterator stop;
	int sz;
public:
	inline bucket() {
		start.ptr = nullptr;
		clear();
	}
	inline ~bucket() {
		clear();
	}
	inline void clear() {
		cup *ptr = start.ptr;
		while (ptr) {
			cup *next = ptr->next;
			bucket_cup_free(ptr);
			ptr = next;
		}
		sz = 0;
		start.ptr = stop.ptr = nullptr;
		start.off = stop.off = 0;
	}
	inline void insert(const particle &part) {
		if (stop.ptr == nullptr) {
			cup *new_cup = bucket_cup_alloc();
			new_cup->next = nullptr;
			new_cup->prev = nullptr;
			start.ptr = new_cup;
			stop.ptr = new_cup;
			stop.off = 0;
		} else if (stop.off == CUP_SIZE) {
			assert(stop.ptr);
			cup *new_cup = bucket_cup_alloc();
			new_cup->next = nullptr;
			new_cup->prev = stop.ptr;
			stop.ptr->next = new_cup;
			stop.ptr = new_cup;
			stop.off = 0;
		}
		stop.ptr->data[stop.off] = part;
		stop.off++;
		sz++;
	}
	inline unsigned size() const {
		return sz;
	}
	inline iterator begin() {
		return start;
	}
	inline iterator end() {
		return stop;
	}
	inline const_iterator cbegin() const {
		return start;
	}
	inline const_iterator cend() const {
		return stop;
	}
	inline iterator remove(iterator iter) {
		assert(start.ptr);
		assert(sz >= 0);
		stop.off--;
		*iter = *stop;
		sz--;
		if (stop.off == 0) {
			assert(stop.ptr);
			auto *free_ptr = stop.ptr;
			if (stop.ptr == start.ptr) {
				start.ptr = stop.ptr = nullptr;
				start.off = stop.off = 0;
			} else {
				stop.ptr = stop.ptr->prev;
				if (stop.ptr->prev) {
					stop.ptr->prev->next = stop.ptr;
				}
				stop.ptr->next = nullptr;
				stop.off = CUP_SIZE;
			}
			if (iter.ptr == free_ptr && iter.off == 0) {
				iter.ptr = stop.ptr;
				iter.off = stop.off;
			}
			bucket_cup_free(free_ptr);
		}
		return iter;
	}
	inline bucket& operator=(const bucket &other) {
		clear();
		for (auto i = other.cbegin(); i != other.cend(); i++) {
			insert(*i);
		}
		return *this;
	}
	inline bucket& operator=(bucket &&other) {
		clear();
		start = other.start;
		stop = other.stop;
		sz = other.sz;
		other.start.ptr = other.stop.ptr = nullptr;
		other.start.off = other.stop.off = 0;
		other.sz = 0;
		return *this;
	}
	inline bucket(const bucket &other) {
		start.ptr = nullptr;
		*this = other;
	}
	inline bucket(bucket &&other) {
		start.ptr = nullptr;
		*this = std::move(other);
	}
	inline particle& front() {
		return *(begin());
	}
	HPX_SERIALIZATION_SPLIT_MEMBER();
	template<class A>
	void load(A &&arc, unsigned) {
		assert(sz==0);
		particle p;
		int this_sz;
		arc & this_sz;
	//	printf( "load %i\n", this_sz);
		for (int i = 0; i < this_sz; i++) {
			arc & p;
			insert(p);
		}
	}
	template<class A>
	void save(A &&arc, unsigned) const {
//		printf( "save %i\n", sz);
		arc & sz;
		for (auto iter = start; iter != stop; iter++) {
			arc & *iter;
		}
	}
};

//struct bucket {
//private:
//	std::vector<particle> v;
//public:
//	inline bucket() {
//	}
//	inline ~bucket() {
//	}
//	inline void clear() {
//		v.clear();
//	}
//	inline void insert(const particle &part) {
//		v.push_back(part);
//	}
//	inline unsigned size() const {
//		return v.size();
//	}
//	inline auto begin() {
//		return v.begin();
//	}
//	inline auto end() {
//		return v.end();
//	}
//	inline auto cbegin() const {
//		return v.cbegin();
//	}
//	inline auto cend() const {
//		return v.cend();
//	}
//	inline auto remove(std::vector<particle>::iterator iter) {
//		*iter = v[v.size() - 1];
//		v.resize(v.size() - 1);
//		if( v.size() == 0 ) {
//			v = decltype(v)();
//		}
//		return iter;
//	}
//	inline bucket& operator=(const bucket &other) {
//		v = other.v;
//		return *this;
//	}
//	inline bucket& operator=(bucket &&other) {
//		v = std::move(other.v);
//		return *this;
//	}
//	inline bucket(const bucket &other) {
//		*this = other;
//	}
//	inline bucket(bucket &&other) {
//		*this = std::move(other);
//	}
//	inline particle& front() {
//		return *(begin());
//	}
//	HPX_SERIALIZATION_SPLIT_MEMBER();
//	template<class A>
//	void load(A &&arc, unsigned) {
//		arc & v;
//	}
//	template<class A>
//	void save(A &&arc, unsigned) const {
//		arc & v;
//	}
//};

#endif /* COSMICTIGER_BUCKET_HPP_ */
