/*
 * timer.hpp
 *
 *  Created on: Sep 18, 2020
 *      Author: dmarce1
 */

#ifndef TIGERGRAV_TIMER_HPP_
#define TIGERGRAV_TIMER_HPP_

#include <string>
#include <chrono>

inline double timer(void) {
	return std::chrono::duration_cast < std::chrono::milliseconds > (std::chrono::system_clock::now().time_since_epoch()).count() / 1000.0;
}



#endif /* TIGERGRAV_TIMER_HPP_ */
