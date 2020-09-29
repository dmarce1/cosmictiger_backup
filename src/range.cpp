/*
 * math.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: dmarce1
 */

#include <cosmictiger/range.hpp>
#include <cosmictiger/util.hpp>

vect<double> range_center(const range &r) {
	vect<double> c;
	for (int dim = 0; dim < NDIM; dim++) {
		c[dim] = (r.min[dim] + r.max[dim]) / 2.0;
	}
	return c;
}

int range_sibling_index(const range &r1, const range &r2) {
	int cnt = 0;
	int rc;
	for (int dim = 0; dim < NDIM, cnt < 2; dim++) {
		if ((r2.min[dim] == r1.max[dim]) || (r2.min[dim] + 1.0 == r1.max[dim])) {
			rc = 2 * dim + 1;
			cnt++;
		} else if ((r1.min[dim] == r2.max[dim]) || (r1.min[dim] + 1.0 == r2.max[dim])) {
			rc = 2 * dim;
			cnt++;
		}
	}
	if (cnt != 1) {
		rc = -1;
	}
	return rc;
}

int range_max_dim(const range &r) {
	int max_dim;
	double max_span = 0.0;
	for (int dim = 0; dim < NDIM; dim++) {
		const auto span = r.max[dim] - r.min[dim];
		if (span > max_span) {
			max_dim = dim;
			max_span = span;
		}
	}
	return max_dim;
}

double range_max_span(const range &r) {
	double s = 0.0;
	for (int dim = 0; dim < NDIM; dim++) {
		s += std::pow(r.max[dim] - r.min[dim], 2);
	}
	return std::sqrt(s);
}

range reflect_range(const range &r_, int dim, double x) {
	range r = r_;
	r.min[dim] = 2 * x - r_.max[dim];
	r.max[dim] = 2 * x - r_.min[dim];
	return r;
}

range scale_range(const range &r_, double s) {
	range r;
	for (int dim = 0; dim < NDIM; dim++) {
		r.min[dim] = r_.min[dim] * s;
		r.max[dim] = r_.max[dim] * s;
	}
	return r;
}

range shift_range(const range &r_, const vect<double> &v) {
	range r;
	for (int dim = 0; dim < NDIM; dim++) {
		r.min[dim] = r_.min[dim] + v[dim];
		r.max[dim] = r_.max[dim] + v[dim];
	}
	return r;
}

bool in_range(const vect<double> &x, const range &r) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (x[dim] < r.min[dim] || x[dim] > r.max[dim]) {
			return false;
		}
	}
	return true;
}

vect<double> range_span(const range &r) {
	vect<double> s;
	for (int dim = 0; dim < NDIM; dim++) {
		s[dim] = r.max[dim] - r.min[dim];
	}
	return s;
}

bool in_range(const range &a, const range &b) {
	return in_range(a.min, b) && in_range(a.max, b);
}

double range_volume(const range &r) {
	double v = 1.0;
	for (int dim = 0; dim < NDIM; dim++) {
		v *= r.max[dim] - r.min[dim];
	}
	return v;
}

bool ranges_intersect(const range &a, const range &b) {
	for (int dim = 0; dim < NDIM; dim++) {
		const auto front = std::max(a.min[dim], b.min[dim]);
		const auto back = std::min(a.max[dim], b.max[dim]);
		if (front > back) {
			return false;
		}
	}
	return true;
}

range null_range() {
	range null;
	for (int dim = 0; dim < NDIM; dim++) {
		null.min[dim] = std::numeric_limits<double>::max();
		null.max[dim] = -std::numeric_limits<double>::max();
	}
	return null;
}

range range_around(const vect<double> &p, double h) {
	range r;
	for (int dim = 0; dim < NDIM; dim++) {
		r.max[dim] = p[dim] + h;
		r.min[dim] = p[dim] - h;
	}
	return r;
}

range range_expand(const range &rs, double h) {
	range rb;
	for (int dim = 0; dim < NDIM; dim++) {
		rb.min[dim] = rs.min[dim] - h;
		rb.max[dim] = rs.max[dim] + h;
	}
	return rb;
}

bool operator==(const range &a, const range &b) {
	for (int dim = 0; dim < NDIM; dim++) {
		if (a.min[dim] != b.min[dim]) {
			return false;
		}
		if (a.max[dim] != b.max[dim]) {
			return false;
		}
	}
	return true;
}

bool operator!=(const range &a, const range &b) {
	return !(a == b);
}

int box_id_to_level(box_id_type id) {
	return msb(id) - 1;
}

range box_id_to_range(box_id_type id) {
	vect<double> x = 0.0;
	vect<double> dx = 1.0;

	while (id != 1) {
//		printf( "%lli\n", id);
		if (id == 0) {
			printf("Logic error line %i file %s\n", __LINE__, __FILE__);
			abort();
		}
		const auto tmp1 = x[0];
		const auto tmp2 = dx[0];
		x[0] = x[2];
		x[2] = x[1];
		x[1] = tmp1;
		dx[0] = dx[2];
		dx[2] = dx[1];
		dx[1] = tmp2;
		if (id & 1) {
			x[0] = 0.5 * x[0] + 0.5;
		} else {
			x[0] = 0.5 * x[0];
		}
		dx[0] *= 0.5;
		id >>= 1;
	}
	range r;
	r.min = x;
	for (int dim = 0; dim < NDIM; dim++) {
		r.max[dim] = r.min[dim] + dx[dim];
	}
	return r;

}

