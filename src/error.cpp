/*
 * error.cpp
 *
 *  Created on: Oct 10, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/error.hpp>
#include <cosmictiger/gravity_cuda.hpp>
#include <unordered_set>

#define SAMPLE_SIZE 100

std::pair<double, double> compute_error(std::vector<output_part> &parts) {
	srand(1234);
	std::pair<double, double> rc;
	std::vector<_4force> f;
	std::vector<vect<double>> x;
	std::vector<vect<double>> y;
	std::vector<int> indices;
	const int sample_size = std::min(SAMPLE_SIZE, (int) parts.size());
	if (sample_size == SAMPLE_SIZE) {
		std::unordered_set<int> i1;
		while (i1.size() < sample_size) {
			i1.insert(rand() % parts.size());
		}
		indices.insert(indices.begin(), i1.begin(), i1.end());
	} else {
		for (int i = 0; i < parts.size(); i++) {
			indices.push_back(i);
		}
	}

	for (int i = 0; i < indices.size(); i++) {
		x.push_back(parts[indices[i]].x);
	}
	f.resize(x.size());
	for (int i = 0; i < parts.size(); i++) {
		y.push_back(parts[i].x);
	}

	gravity_ewald_direct(f, x, y);

	double err = 0.0;
	vect<double> gtot = 0.0;
	double gnorm = 0.0;
	for (int i = 0; i < parts.size(); i++) {
		gtot += parts[i].g;
		gnorm += abs(parts[i].g);
	}
	for (int i = 0; i < f.size(); i++) {
		auto dg = parts[indices[i]].g - f[i].g;
	//	parts[indices[i]].g =f[i].g;
		err += abs(dg)/(abs(f[i].g));
	//	printf( "%e %e\n", parts[indices[i]].g[0], f[i].g[0]);
	}
	err /= f.size();
	gtot /= gnorm;
	printf( "%.14e %.14e %.14e %.14e \n", err, gtot[0], gtot[1], gtot[2]);
	return rc;
}
