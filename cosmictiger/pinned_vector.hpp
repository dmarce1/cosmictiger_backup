/*
 * pinned_vector.hpp
 *
 *  Created on: Sep 18, 2020
 *      Author: dmarce1
 */

#ifndef TIGERGRAV_PINNED_VECTOR_HPP_
#define TIGERGRAV_PINNED_VECTOR_HPP_


#include <cosmictiger/cuda_check.hpp>

#include <cuda_runtime.h>

template<class T>
class pinned_vector {
	T *ptr;
	std::size_t sz;
	std::size_t cap;
	void free() {
		if (ptr) {
			cudaFreeHost(ptr);
		}
	}
	void allocate() {
		CUDA_CHECK(cudaMallocHost((void**) &ptr, sizeof(T) * cap));
	}
public:
	pinned_vector() {
		cap = 0;
		sz = 0;
		ptr = nullptr;
	}
	pinned_vector(const pinned_vector&) = delete;
	pinned_vector& operator=(const pinned_vector&) = delete;
	pinned_vector(std::size_t this_sz) {
		cap = this_sz;
		sz = this_sz;
		allocate();
	}
	~pinned_vector() {
		free();
	}
	pinned_vector& operator=(pinned_vector &&other) {
		free();
		cap = other.cap;
		sz = other.sz;
		ptr = other.ptr;
		other.cap = 0;
		other.sz = 0;
		other.ptr = nullptr;
		return *this;
	}
	pinned_vector(pinned_vector &&other) {
		ptr = nullptr;
		(*this) = std::move(other);
	}
	void resize(std::size_t new_sz) {
		if (new_sz > cap) {
			T *new_ptr;
			CUDA_CHECK(cudaMallocHost((void**) &new_ptr, sizeof(T) * new_sz));
			if (sz) {
				std::memcpy(new_ptr, ptr, sz * sizeof(T));
			}
			free();
			ptr = new_ptr;
			cap = new_sz;
			sz = new_sz;
		} else {
			sz = new_sz;
		}
	}
	void reserve(std::size_t new_cap) {
		if (cap < new_cap) {
			auto old_sz = sz;
			resize(new_cap);
			sz = old_sz;
		}
	}
	void push_back(T &&v) {
		if (sz + 1 >= cap) {
			reserve(std::max(2 * cap, (std::size_t) 1));
		}
		ptr[sz] = std::move(v);
		sz++;
	}
	void push_back(const T &v) {
		if (sz + 1 >= cap) {
			reserve(std::max(2 * cap, (std::size_t) 1));
		}
		ptr[sz] = v;
		sz++;
	}
	T* data() {
		return ptr;
	}
	const T* data() const {
		return ptr;
	}
	std::size_t size() const {
		return sz;
	}
	T operator[](int i) const {
		return ptr[i];
	}
	T& operator[](int i) {
		return ptr[i];
	}
	T* begin() {
		return ptr;
	}
	T* end() {
		return ptr + sz;
	}
//	template<class I>
//	void insert(T *start, I begin, I end) {
//		printf("Insert begin\n");
//		const std::size_t rsz = (start - ptr) + (end - begin);
//		if (rsz > cap) {
//			reserve(rsz);
//		}
//		for (auto i = begin; i != end; i++) {
//			*(start + std::distance(i, begin)) = *i;
//		}
//		sz = std::max(sz, rsz);
//		printf("Insert end\n");
//	}
};



#endif /* TIGERGRAV_PINNED_VECTOR_HPP_ */
