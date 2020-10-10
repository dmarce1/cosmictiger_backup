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
	part_pos x;
	vect<float> v;
	struct {
		std::uint64_t rung :7;
		std::uint64_t out :1;
		std::uint64_t group :56;
	};
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & x;
		arc & v;
		std::uint8_t tmp1;
		std::uint64_t tmp2;
		tmp1 = rung;
		arc & tmp1;
		rung = tmp1;
		tmp1 = out;
		arc & tmp1;
		out = tmp1;
		tmp2 = group;
		arc & tmp2;
		group = tmp2;
	}
};

#define DEFAULT_GROUP 0

#endif /* COSMICTIGER_PARTICLE_HPP_ */
