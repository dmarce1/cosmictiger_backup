/*
 * gravity.hpp
 *
 *  Created on: Oct 4, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_GRAVITY_HPP_
#define COSMICTIGER_GRAVITY_HPP_

#include <cosmictiger/expansion.hpp>
#include <cosmictiger/vect.hpp>

struct _4force {
	float phi;
	vect<float> g;
};

#define SELF_PHI 2.8

std::uint64_t gravity_CP_direct(expansion<float> &L, const vect<position> &x, std::vector<vect<position>> y, bool do_phi);
std::uint64_t gravity_CC_direct(expansion<float> &L, const vect<position> &x, std::vector<const multi_src*> &y, bool do_phi);
std::uint64_t gravity_CC_ewald(expansion<float> &L, const vect<position> &x, std::vector<const multi_src*> &y, bool do_phi);


#endif /* COSMICTIGER_GRAVITY_HPP_ */
