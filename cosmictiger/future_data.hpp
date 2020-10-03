/*
 * future_data.hpp
 *
 *  Created on: Oct 3, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_FUTURE_DATA_HPP_
#define COSMICTIGER_FUTURE_DATA_HPP_

#include <cosmictiger/hpx.hpp>

template<class T>
class future_data {
	hpx::future<T> fut;
	T data;
public:
	void set(const T &d) {
		data = d;
	}
	void set(hpx::future<T> &&f) {
		fut = std::move(f);
	}
	T get() {
		if (fut.valid()) {
			return fut.get();
		} else {
			return data;
		}
	}
	bool is_ready() const {
		if (fut.valid()) {
			return fut.is_ready();
		} else {
			return true;
		}
	}
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & data;
		arc & fut;
	}
};

template<class T>
class shared_future_data {
	hpx::shared_future<T> fut;
	T data;
public:
	void set(const T &d) {
		data = d;
	}
	void set(hpx::future<T> &&f) {
		fut = std::move(f);
	}
	T& get() {
		if (fut.valid()) {
			data = fut.get();
		}
		return data;
	}
	bool is_ready() const {
		if (fut.valid()) {
			return fut.is_ready();
		} else {
			return true;
		}
	}
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & data;
		arc & fut;
	}
};

#endif /* COSMICTIGER_FUTURE_DATA_HPP_ */
