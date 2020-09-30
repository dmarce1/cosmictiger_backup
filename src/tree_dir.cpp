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
	level = msb(hpx_localities().size()) - 1;
	nx = 1 << level;
	ny = 1 << (level - 1);
	nz = 1 << (level - 2);
}

void tree_dir::add_tree_client(const tree_client &client, const range &box) {
	nodes[index(range_center(box))] = client;
}

tree_dir& tree_dir::merge(const tree_dir &other) {
	for (auto i : other.nodes) {
		nodes[i.first] = i.second;
	}
	return *this;
}

int tree_dir::index(const vect<double> &x) const {
	return index((int) x[0] * nx, (int) x[1] * ny, (int) x[2] * nz);
}

int tree_dir::index(int x, int y, int z) const {
	return x * ny * nz + y * nz + z;
}

hpx::future<void> tree_dir::find_home(int stack_cnt, bucket &&parts) {
	using map_type = std::unordered_map<int, bucket>;
	map_type buckets;
	while (parts.size()) {
		const auto &p = parts.front();
		const auto i = index(pos_to_double(p.x));
		buckets[i].insert(p);
		parts.remove(parts.begin());
	}
	std::vector<hpx::future<int>> futs;
	for (auto &i : buckets) {
		auto func = [this, stack_cnt](int id, bucket&& parts) {
			nodes[id].find_home_child(stack_cnt + 1, std::move(parts));
			return 0;
		};
		futs.push_back(hpx::async(std::move(func), i.first, std::move(i.second)));
	}
	return hpx::when_all(futs.begin(), futs.end());
}

