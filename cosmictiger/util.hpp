/*
 * util.hpp
 *
 *  Created on: Sep 29, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_UTIL_HPP_
#define COSMICTIGER_UTIL_HPP_

#include <cstdint>

inline int msb(int i) {
	return 63 - __builtin_clzll(i);
}

inline int bits_to_level(int i) {
	return msb(i - 1) + 1;
}
#endif /* COSMICTIGER_UTIL_HPP_ */
