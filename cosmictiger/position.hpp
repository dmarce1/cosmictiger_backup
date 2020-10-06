#pragma once
#include <cosmictiger/defs.hpp>
#include <cosmictiger/vect.hpp>

#include <cstdint>
#include <limits>

static constexpr double POS_MAX = double(std::numeric_limits < std::uint32_t > ::max()) + 1.0;
static constexpr double POS_INV = 1.0 / POS_MAX;

class position {
	std::int32_t i;
public:
	/**/DEFAULT_CLASS_MEMBERS(position);
	position(double r) {
		while (r >= 1.0) {
			r -= 1.0;
		}
		while (r < 0.0) {
			r += 1.0;
		}
		i = (r - 0.5) * POS_MAX;
//		printf( "%e %i %e\n", r, i, POS_INV);
	}
	operator double() const {
		return double(i) * POS_INV + 0.5;
	}
	operator int() const {
		return i;
	}
	position operator-(const position &other) const {
		position rc;
		rc.i = i - other.i;
		return rc;
	}
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & i;
	}
};

inline double pos_to_double(const position &p1) {
	return double(p1);
}

inline vect<double> pos_to_double(const vect<position> &p1) {
	vect<double> p2;
	for (int dim = 0; dim < NDIM; dim++) {
		p2[dim] = double(p1[dim]);
	}
	return p2;
}

inline vect<double> double_to_pos(const vect<double> &p1) {
	vect<position> p2;
	for (int dim = 0; dim < NDIM; dim++) {
		p2[dim] = p1[dim];
	}
	return p2;
}

using part_pos = vect<position>;

