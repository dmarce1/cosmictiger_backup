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

#define NSTACK 8

static std::stack<cup*> stack[NSTACK];

#define CHUNK_SIZE (1024*1024)

static std::atomic<int> lock[NSTACK] = {0};

static std::atomic<int> next(0);

//static std::atomic<int> cnt(0);

cup* bucket_cup_alloc() {
	next++;
	int si = next % NSTACK;
	while (lock[si]++ != 0) {
		lock[si]--;
	}
	if (stack[si].empty()) {
//		printf( "Allocating\n");
		auto *base = (cup*) malloc(sizeof(cup) * CHUNK_SIZE);
		for (int i = 0; i < CHUNK_SIZE; i++) {
			stack[si].push(base + i);
		}
	}
	auto *ptr = stack[si].top();
	stack[si].pop();
	lock[si]--;
//	cnt++;
//	if( cnt % 100000 == 0 ) {
//		printf( "%i\n", (int) cnt);
//	}
	return ptr;
}

void bucket_cup_free(cup *ptr) {
	next++;
	int si = next % NSTACK;
	while (lock[si]++ != 0) {
		lock[si]--;
	}
	stack[si].push(ptr);
	lock[si]--;
//	cnt--;
//	if( cnt % 100000 == 0 ) {
//		printf( "%i\n", (int) cnt);
//	}
}
