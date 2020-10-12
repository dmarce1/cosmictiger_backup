#pragma once

#include <cosmictiger/check.hpp>

std::vector<std::shared_ptr<std::vector<part_pos>>> get_positions(const std::vector<tree_ptr> &ids);
void pos_cache_cleanup();
