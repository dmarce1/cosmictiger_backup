#include <cosmictiger/defs.hpp>
#include <cosmictiger/tree.hpp>

#include <atomic>
#include <stack>
#include <thread>

#define MAX_STACK 16

HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<tree>, tree);

static std::stack<tree_mems*> stack;
static mutex_type mtx;

#define CHUNK_SIZE (512*1024)

tree_mems* tree_mems_alloc() {
	std::lock_guard<mutex_type> lock(mtx);
	if (stack.empty()) {
		tree_mems *base = (tree_mems*) malloc(sizeof(tree_mems) * CHUNK_SIZE);
		for (int i = 0; i < CHUNK_SIZE; i++) {
			stack.push(base + i);
		}
	}
	tree_mems *ptr = stack.top();
	stack.pop();
	return ptr;
}

void tree_mems_free(tree_mems *ptr) {
	std::lock_guard<mutex_type> lock(mtx);
	stack.push(ptr);
}

tree::tree() {
	tptr = tree_mems_alloc();
}

tree::tree(box_id_type id) {
	tptr = tree_mems_alloc();
	tptr->boxid = id;
	tptr->leaf = true;
}

tree::~tree() {
	tree_mems_free(tptr);
}

void tree::add_parts(std::vector<particle> &&parts, int stack_cnt) {
	if (parts.size() <= opts.bucket_size) {
		for (const auto &part : parts) {
			tptr->parts.push(part);
		}
	} else {
		create_children();
	}
}

void tree::create_children() {
	auto futl = hpx::new_ < tree > (hpx::find_here(), (tptr->boxid << 1) + 0);
	auto futr = hpx::new_ < tree > (hpx::find_here(), (tptr->boxid << 1) + 1);
	tptr->leaf = false;
	tptr->children[0] = futl.get();
	tptr->children[1] = futr.get();
}

tree_client::tree_client(hpx::id_type &&myid) {
	id = std::move(myid);
}

static bool thread_me(int stack_cnt, hpx::id_type gid) {
	static std::atomic<int> thread_cnt(0);
	static const int max_threads = 4 * std::thread::hardware_concurrency();
	bool threadme;
	if (thread_cnt++ < max_threads) {
		threadme = true;
	} else {
		thread_cnt--;
		threadme = false;
	}
	if (!threadme) {
		if (gid != hpx::find_here() || stack_cnt >= MAX_STACK) {
			thread_cnt++;
			threadme = true;
		}
	}
	return threadme;
}

hpx::future<void> tree_client::add_parts(std::vector<particle> &&parts, int stack_cnt) {
	stack_cnt++;
	if (thread_me(stack_cnt, hpx::get_colocation_id(id).get())) {
		return hpx::async < tree::add_parts_action > (id, std::move(parts), stack_cnt);
	} else {
		tree::add_parts_action()(id, std::move(parts), stack_cnt);
		return hpx::make_ready_future();
	}
}

