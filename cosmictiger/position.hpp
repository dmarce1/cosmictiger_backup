#pragma once
#include <cosmictiger/defs.hpp>
#include <cosmictiger/vect.hpp>

#include <cstdint>
#include <limits>

class position {
	static constexpr double pos_max = double(std::numeric_limits < std::int32_t > ::max()) + 1.0;
	static constexpr double pos_inv = 1.0 / pos_max;
	std::int32_t i;
public:DEFAULT_CLASS_MEMBERS(position)
	;
	position(double r) {
		while (r >= 1.0) {
			r -= 1.0;
		}
		while (r < 0.0) {
			r += 1.0;
		}
		i = (r - 0.5) * pos_max;
	}
	operator double() const {
		return (double(i) + 0.5) * pos_inv + 0.5;
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
