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

static std::stack<cup*> main_stack;
static mutex_type mtx;


class cup_stack {
	std::stack<cup*> stack;
public:
	bool empty() const {
		return stack.empty();
	}
	void pop() {
		stack.pop();
	}
	void push(cup* ptr) {
		stack.push(ptr);
	}
	cup* top() {
		return stack.top();
	}
	std::size_t size() const {
		return stack.size();
	}
	~cup_stack() {
		std::lock_guard<mutex_type> lock(mtx);
		while( !empty()) {
			main_stack.push(top());
			pop();
		}
	}
};

static thread_local cup_stack stack;

#define CHUNK_SIZE (1024)

cup* bucket_cup_alloc() {
	if (stack.empty()) {
		std::lock_guard<mutex_type> lock(mtx);
		if( main_stack.size() < CHUNK_SIZE) {
			auto *base = (cup*) malloc(sizeof(cup) * CHUNK_SIZE);
			for (int i = 0; i < CHUNK_SIZE; i++) {
				main_stack.push(base + i);
			}
		}
		for( int i = 0; i < CHUNK_SIZE; i++) {
			stack.push(main_stack.top());
			main_stack.pop();
		}
	}
	auto *ptr = stack.top();
	stack.pop();
	return ptr;
}

void bucket_cup_free(cup *ptr) {
	stack.push(ptr);
	if (stack.size() >= 2*CHUNK_SIZE) {
		std::lock_guard<mutex_type> lock(mtx);
		for( int i = 0; i < CHUNK_SIZE; i++) {
			main_stack.push(stack.top());
			stack.pop();
		}
	}
}
