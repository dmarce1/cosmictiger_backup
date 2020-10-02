#pragma once


#include <string>

void fileio_init_read(const std::string basename);
std::uint64_t fileio_insert_parts(std::uint64_t chunk_size);


