#pragma once

#include <cosmictiger/gravity.hpp>
#include <cosmictiger/multipole.hpp>
#include <cosmictiger/tree_client.hpp>

#include <functional>
#include <memory>
#include <vector>

std::uint64_t gravity_queue_genid();
void gravity_queue_checkin(std::uint64_t);
void gravity_queue_retire_futures();

void gravity_queue_add_work(std::uint64_t, std::shared_ptr<std::vector<_4force>>, std::shared_ptr<std::vector<part_pos>>, std::vector<tree_client>&&,
		std::vector<multi_src*>&&, std::function<void(void)>&&);
