/*
 * particle.hpp
 *
 *  Created on: Sep 25, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_PARTICLE_HPP_
#define COSMICTIGER_PARTICLE_HPP_

#include <cosmictiger/position.hpp>
#include <cosmictiger/vect.hpp>

struct particle {
	vect<position> x;
	vect<float> v;
	struct {
		std::uint64_t rung :7;
		std::uint64_t out :1;
		std::uint64_t group :56;
	};
};

#endif /* COSMICTIGER_PARTICLE_HPP_ */
