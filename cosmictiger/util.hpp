/*
 * util.hpp
 *
 *  Created on: Sep 29, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_UTIL_HPP_
#define COSMICTIGER_UTIL_HPP_

#include <cstdint>

inline int msb(std::uint64_t i) {
	return 64 - __builtin_clz(i);
}

#endif /* COSMICTIGER_UTIL_HPP_ */
