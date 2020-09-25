/*
 * hpx.cpp
 *
 *  Created on: Sep 25, 2020
 *      Author: dmarce1
 */

#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>

static mutex_type mtx;
static std::vector<hpx::id_type> localities;
static bool initialized = false;

const std::vector<hpx::id_type>& hpx_localities() {
	if (!initialized) {
		std::lock_guard<mutex_type> lock(mtx);
		if (!initialized) {
			localities = hpx::find_all_localities();
			initialized = true;
		}
	}
	return localities;
}

