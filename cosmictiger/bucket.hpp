/*
 * bucket.hpp
 *
 *  Created on: Sep 25, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_BUCKET_HPP_
#define COSMICTIGER_BUCKET_HPP_

#include <cosmictiger/defs.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/particle.hpp>

#include <cassert>

class bucket;

particle* bucket_alloc();
void bucket_free(particle *ptr);

class bucket {
	static const options opts;
	particle *parts;
	unsigned sz;
public:
	inline bucket() {
		sz = 0;
	}
	inline ~bucket() {
		if (parts) {
			bucket_free(parts);
		}
	}
	inline particle* data() {
		return parts;
	}
	inline particle operator[](unsigned i) const {
		assert(i >= 0);
		assert(i < sz);
		return parts[i];
	}
	inline particle& operator[](unsigned i) {
		assert(i >= 0);
		assert(i < sz);
		return parts[i];
	}
	inline void resize(unsigned new_size) {
		assert(new_size < opts.bucket_size);
		sz = new_size;
		if (sz > 0 && !parts) {
			parts = bucket_alloc();
		} else if (sz == 0 && parts) {
			bucket_free(parts);
		}
	}
	inline void push(const particle &part) {
		if (!parts) {
			parts = bucket_alloc();
		}
		parts[sz] = part;
		sz++;
		assert(sz < opts.bucket_size);
	}
	inline unsigned size() const {
		return sz;
	}
	inline particle* begin() {
		return parts;
	}
	inline particle* end() {
		return parts + sz;
	}
	inline void remove(unsigned i) {
		assert(i < sz);
		assert(sz > 0);
		parts[i] = parts[sz - 1];
		sz--;
		if( sz == 0 ) {
			bucket_free(parts);
		}
	}
	inline void remove(particle *iter) {
		remove(iter - parts);
	}
	inline bucket& operator=(const bucket &other) {
		resize(other.sz);
		for (unsigned i = 0; i < sz; i++) {
			parts[i] = other.parts[i];
		}
		return *this;
	}
	inline bucket& operator=(bucket &&other) {
		sz = other.sz;
		parts = other.parts;
		other.sz = 0;
		if( other.parts) {
			bucket_free(other.parts);
		}
		return *this;
	}
	inline bucket(const bucket &other) {
		*this = other;
	}
	inline bucket(bucket &&other) {
		*this = std::move(other);
	}

	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & sz;
		if (sz > 0 && !parts) {
			parts = bucket_alloc();
		} else if (sz == 0 && parts) {
			bucket_free(parts);
		}
		for (int i = 0; i < sz; i++) {
			arc & parts[i];
		}
	}
};

#endif /* COSMICTIGER_BUCKET_HPP_ */
