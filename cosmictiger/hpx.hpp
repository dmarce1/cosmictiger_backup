#pragma once


#ifdef HPX_LITE
#include <hpx/hpx_lite.hpp>
#else
#include <hpx/hpx_init.hpp>
#include <hpx/async.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/plain_actions.hpp>
#include <hpx/runtime/get_colocation_id.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/include/lcos.hpp>
#endif

using mutex_type = hpx::lcos::local::spinlock;

const std::vector<hpx::id_type>& hpx_localities();
