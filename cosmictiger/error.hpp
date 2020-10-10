/*
 * error.hpp
 *
 *  Created on: Oct 10, 2020
 *      Author: dmarce1
 */

#ifndef COSMICTIGER_ERROR_HPP_
#define COSMICTIGER_ERROR_HPP_

#include <cosmictiger/output.hpp>

std::pair<double, double> compute_error(const std::vector<output_part>&);

#endif /* COSMICTIGER_ERROR_HPP_ */
