/*
 * math.hpp
 *
 *  Created on: Nov 29, 2019
 *      Author: dmarce1
 */

#ifndef RANGE_HPP_
#define RANGE_HPP_

#include <cosmictiger/vect.hpp>

#include <limits>

struct range {
	vect<double> min;
	vect<double> max;
	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & min;
		arc & max;
	}
};

using box_id_type = std::uint64_t;

range box_id_to_range(box_id_type id);
int box_id_to_level(box_id_type id);

int range_max_dim(const range& r );
int range_sibling_index(const range &r1, const range &r2);
double range_max_span(const range &r);
range reflect_range(const range&, int dim, double axis);
vect<double> range_center(const range &r);
range shift_range(const range &r, const vect<double>&);
range scale_range(const range&, double);
vect<double> range_span(const range&);
bool in_range(const vect<double>&, const range&);
bool in_range(const range&, const range&);
bool ranges_intersect(const range&, const range&);
range range_around(const vect<double>&, double);
range range_expand(const range&, double);
double range_volume(const range&);
range null_range();
bool operator==(const range&, const range&);
bool operator!=(const range&, const range&);


inline bool in_range(const vect<double> &x, const range &r) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (x[dim] < r.min[dim] || x[dim] > r.max[dim]) {
			return false;
		}
	}
	return true;
}
#endif /* MATH_HPP_ */

