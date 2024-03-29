#pragma once

#include <cosmictiger/options.hpp>

#include <limits>

using rung_type = std::int8_t;
using time_type = std::uint32_t;

#define RUNG_MAX 30

inline double time_to_double(time_type t) {
	const auto t_max = opts.t_max;
	const auto imax = time_type(1) << time_type(31);
	return ((double) t / (double) imax) * t_max;
}

inline time_type inc(time_type t, rung_type max_rung) {
	t += (time_type(1) << time_type(31 - max_rung));
	return t;
}

inline rung_type min_rung(time_type t) {
	rung_type min_rung = 31;
	while (((t & 1) == 0) && (min_rung != 0)) {
		min_rung--;
		t >>= 1;
	}
	return min_rung;
}

inline double rung_to_dt(std::int8_t rung) {
	const auto t_max = opts.t_max;
	return t_max / (1 << rung);
}

inline rung_type dt_to_rung(double dt) {
	int rung = 0;
	while (rung_to_dt(rung) > dt) {
//		printf( "%i %e\n", rung, dt);
		rung++;
		if (rung == std::numeric_limits < rung_type > ::max()) {
			printf("logic error %s %i\n", __FILE__, __LINE__);
			abort();
		}
	}
	return rung;
}

