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
#include <vector>

static std::stack<cup*, std::vector<cup*>> main_stack;
static std::atomic<int> lock(0);

class cup_stack {
	std::stack<cup*, std::vector<cup*>> stack;
public:
	bool empty() const {
		return stack.empty();
	}
	void pop() {
		stack.pop();
	}
	void push(cup *ptr) {
		stack.push(ptr);
	}
	cup* top() {
		return stack.top();
	}
	std::size_t size() const {
		return stack.size();
	}
	~cup_stack() {
		while (lock++ != 0) {
			lock--;
		}
		while (!empty()) {
			main_stack.push(top());
			pop();
		}
		lock--;
	}
};

static thread_local cup_stack stack;

#define CHUNK_SIZE (1024)

cup* bucket_cup_alloc() {
	if (stack.empty()) {
		while (lock++ != 0) {
			lock--;
		}
		if (main_stack.size() < CHUNK_SIZE) {
			auto *base = (cup*) malloc(sizeof(cup) * CHUNK_SIZE);
			for (int i = 0; i < CHUNK_SIZE; i++) {
				main_stack.push(base + i);
			}
		}
		for (int i = 0; i < CHUNK_SIZE; i++) {
			stack.push(main_stack.top());
			main_stack.pop();
		}
		lock--;
	}
	auto *ptr = stack.top();
	stack.pop();
	return ptr;
}

void bucket_cup_free(cup *ptr) {
	stack.push(ptr);
	if (stack.size() >= 2 * CHUNK_SIZE) {
		while (lock++ != 0) {
			lock--;
		}
		for (int i = 0; i < CHUNK_SIZE; i++) {
			main_stack.push(stack.top());
			stack.pop();
		}
		lock--;
	}
}
