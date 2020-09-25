#pragma once
#include <cosmictiger/defs.hpp>

#include <cstdint>
#include <limits>



class position {
	static constexpr double pos_max = double(std::numeric_limits<std::int32_t>::max()) + 1.0;
	static constexpr double pos_inv = 1.0 / pos_max;
	std::int32_t i;
public:
	DEFAULT_CLASS_MEMBERS(position);
	position(double r) {
		i = (r - 0.5) * pos_max;
	}
	operator double() const {
		return (double(i) + 0.5) * pos_inv + 0.5;
	}
	position operator-(const position& other) const {
		position rc;
		rc.i = i - other.i;
		return rc;
	}
};
