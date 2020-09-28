/*
 * bucket.cpp
 *
 *  Created on: Sep 25, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/defs.hpp>
#include <cosmictiger/bucket.hpp>
#include <cosmictiger/hpx.hpp>


#include <stack>

static std::stack<cup*> stack;
static mutex_type mtx;

#define CHUNK_SIZE (1024*1024)

cup* bucket_cup_alloc() {
	std::lock_guard < mutex_type > lock(mtx);
	if (stack.empty()) {
		auto *base = (cup*) malloc(sizeof(cup) * CHUNK_SIZE);
		for (int i = 0; i < CHUNK_SIZE; i++) {
			stack.push(base + i);
		}
	}
	auto *ptr = stack.top();
	stack.pop();
	return ptr;
}

void bucket_cup_free(cup *ptr) {
	std::lock_guard < mutex_type > lock(mtx);
	stack.push(ptr);
}
