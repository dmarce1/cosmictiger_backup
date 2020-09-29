/*
 * hpx.cpp
 *
 *  Created on: Sep 25, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>

static std::atomic<int> lock;
static std::vector<hpx::id_type> localities;
static bool initialized = false;

const std::vector<hpx::id_type>& hpx_localities() {
	if (!initialized) {
		while (lock++ != 0) {
			lock--;
		}
		if (!initialized) {
			localities = hpx::find_all_localities();
			initialized = true;
		}
		lock--;
	}
	return localities;
}

