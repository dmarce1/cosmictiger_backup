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

static std::stack<particle*> stack;
static mutex_type mtx;
static bool initialized = false;

#define CHUNK_SIZE (512*1024)

particle* bucket_alloc() {
	std::lock_guard < mutex_type > lock(mtx);
	if (stack.empty()) {
		particle *base = (particle*) malloc(sizeof(particle) * CHUNK_SIZE * opts.bucket_size);
		for (int i = 0; i < CHUNK_SIZE; i++) {
			stack.push(base + i);
		}
	}
	particle *ptr = stack.top();
	stack.pop();
	return ptr;
}

void bucket_free(particle *ptr) {
	std::lock_guard < mutex_type > lock(mtx);
	stack.push(ptr);
}
