#include <cosmictiger/check.hpp>
#include <cosmictiger/hpx.hpp>

#include <stack>

std::stack<void*, std::vector<void*>> stack;
mutex_type mtx;

void* check_allocate(std::size_t size) {
	auto* ptr = malloc(size);
	std::lock_guard<mutex_type> lock(mtx);
	return ptr;
}

void check_cleanup() {
	while (stack.size()) {
		free(stack.top());
		stack.pop();
	}
}

