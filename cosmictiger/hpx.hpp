#pragma once


#ifdef HPX_LITE
#include <hpx/hpx_lite.hpp>
#else
#include <hpx/hpx_init.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/plain_actions.hpp>
#endif

using mutex_type = hpx::lcos::local::spinlock;
