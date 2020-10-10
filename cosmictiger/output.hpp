#pragma once

#include <cosmictiger/gravity.hpp>
#include <cosmictiger/particle.hpp>
#include <cosmictiger/hpx.hpp>

#include <string>

struct output_part {
	vect<double> x;
	vect<float> v;
	vect<float> g;
	float phi;
	int rung;
	std::uint64_t group_id;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & x;
		arc & v;
		arc & g;
		arc & phi;
		arc & rung;
	}
};

void output_add_particle(const particle&, const _4force&);
void output_to_file(const std::string&);
