/*
 * error.cpp
 *
 *  Created on: Oct 10, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/error.hpp>
#include <cosmictiger/gravity_cuda.hpp>

std::pair<double, double> compute_error(const std::vector<output_part> &parts) {
	std::pair<double, double> rc;
	std::vector<_4force> f;
	std::vector<vect<double>> x;
	std::vector<vect<double>> y;
	std::vector<int> indices;
	for (int i = 0; i < parts.size(); i++) {
		x.push_back(parts[i].x);
	}
	f.resize(x.size());
	for (int i = 0; i < parts.size(); i++) {
		y.push_back(parts[i].x);
	}

	gravity_ewald_direct(f, x, y);

	for (int i = 0; i < f.size(); i++) {
//		parts[i].phi = f[i].phi;
//		parts[i].g = f[i].g;
		printf("%e %e\n", f[i].phi, parts[i].phi);
	}

	return rc;
}
