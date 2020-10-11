/*
 * error.cpp
 *
 *  Created on: Oct 10, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/timer.hpp>
#include <cosmictiger/error.hpp>
#include <cosmictiger/gravity_cuda.hpp>
#include <unordered_set>

#define SAMPLE_SIZE 500

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
	auto ts = timer();
	gravity_ewald_direct(f, x, y);
	printf("Direct solve took %e seconds\n", timer() - ts);
	double err = 0.0;
	double err2 = 0.0;
	vect<double> gtot = 0.0;
	double gnorm = 0.0;
	for (int i = 0; i < parts.size(); i++) {
		gtot += parts[i].g;
		gnorm += abs(parts[i].g);
	}
	std::vector<double> errs;
	for (int i = 0; i < f.size(); i++) {
		auto dg = parts[indices[i]].g - f[i].g;
		const auto this_err = abs(dg) / (abs(f[i].g));
		err += this_err;
		err2 += this_err * this_err;
		errs.push_back(this_err);
		//	printf( "%e %e\n", parts[indices[i]].g[0], f[i].g[0]);
	}
	err /= f.size();
	err2 /= f.size();
	err2 = std::sqrt(err2);
	gtot /= gnorm;
	std::sort(errs.begin(), errs.end());
	const auto err99 = errs[SAMPLE_SIZE * 99 / 100];
	printf("Relative Error : %e\n", err);
	printf("RMS      Error : %e\n", err2);
	printf("99th%%    Error : %e\n", err99);
	printf("gtot = %.14e %.14e %.14e \n",  gtot[0], gtot[1], gtot[2]);
	return rc;
}
