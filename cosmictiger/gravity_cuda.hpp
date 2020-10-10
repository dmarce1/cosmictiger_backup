/*
 * gravity_cuda.hpp
 *
 *  Created on: Sep 1, 2020
 *      Author: dmarce1
 */

#ifndef TIGERGRAV_GRAVITY_CUDA_HPP_
#define TIGERGRAV_GRAVITY_CUDA_HPP_

#include <cosmictiger/cuda_check.hpp>
#include <cosmictiger/gravity.hpp>
#include <cosmictiger/particle.hpp>
#include <cosmictiger/pinned_vector.hpp>
#include <cosmictiger/expansion.hpp>


#define EWALD_NFOUR 80
#define EWALD_NREAL 171

double cuda_reset_flop();

#include <cosmictiger/gravity.hpp>

#include <memory>

#define SYNCRATE 10

struct cuda_ewald_const {
	vect<float> four_indices[EWALD_NFOUR];
	vect<float> real_indices[EWALD_NREAL];
	expansion<float> periodic_parts[EWALD_NFOUR];
	expansion<float> exp_factors;
};

struct cuda_work_unit {
	std::vector<std::pair<int, int>> yiters;
	std::vector<const multi_src*> z;
	std::shared_ptr<std::vector<vect<position>>> xptr;
	std::shared_ptr<std::vector<_4force>> fptr;
};


template<class T>
class pinned_vector;

void gravity_PP_direct_cuda(std::vector<cuda_work_unit>&&, const pinned_vector<vect<position>>&, bool do_phi);
void gravity_CC_ewald_cuda(expansion<float> &L, const vect<position> &x, std::vector<const multi_src*> &y, bool do_phi);
void gravity_ewald_direct(std::vector<_4force> &f, const std::vector<vect<double>> x, const std::vector<vect<double>> &y);

#endif /* TIGERGRAV_GRAVITY_CUDA_HPP_ */
