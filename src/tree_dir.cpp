/*
 * tree_top.cpp
 *
 *  Created on: Sep 29, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/hpx.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/tree_dir.hpp>
#include <cosmictiger/util.hpp>

tree_dir::tree_dir() {
	level = bits_to_level(hpx_localities().size());
	nx = 1 << ((level + 2) / NDIM);
	ny = 1 << ((level + 1) / NDIM);
	nz = 1 << ((level + 0) / NDIM);
}

void tree_dir::add_tree_client(const tree_client &client, const range &box) {
	const int id = index(range_center(box));
	nodes[id] = client;
	assert((nodes[id] != tree_client()));
}

tree_dir& tree_dir::merge(const tree_dir &other) {
	for (auto i : other.nodes) {
		nodes[i.first] = i.second;
	}
	for (auto i : nodes) {
		assert(i.second != tree_client());
	}
	return *this;
}

int tree_dir::index(const vect<double> &x) const {
	return index((int) (x[0] * nx), (int) (x[1] * ny), (int) (x[2] * nz));
}

int tree_dir::index(int x, int y, int z) const {
	return x * ny * nz + y * nz + z;
}

hpx::future<void> tree_dir::find_home(bucket &&parts) {
	using map_type = std::unordered_map<int, bucket>;
	map_type buckets;
	while (parts.size()) {
		const auto &p = parts.front();
		const auto x = pos_to_double(p.x);
		const auto i = index(x);
//		printf("Sending %e %e %e to %i\n", x[0], x[1], x[2], i);
		buckets[i].insert(p);
		parts.remove(parts.begin());
	}
	std::vector<hpx::future<int>> futs;
	for (auto &i : buckets) {
		auto func = [this](int id, bucket &&parts) {
			assert((nodes[id] != tree_client()));
			nodes[id].place_parts(std::move(parts));
			return 0;
		};
		futs.push_back(hpx::async(std::move(func), i.first, std::move(i.second)));
	}
	std::lock_guard<mutex_type> lock(mtx);
	return hpx::when_all(futs.begin(), futs.end()).then([](hpx::future<std::vector<hpx::future<int>>> fut) {
		auto futs = fut.get();
		hpx::wait_all(futs.begin(), futs.end());
	});
}

