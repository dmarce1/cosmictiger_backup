# Copyright (c) 2015 Thomas Heller
# Copyright (c) 2015 Dominic Marcello
# Copyright (c) 2018-2019 Parsa Amini
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.12.4)
cmake_policy(VERSION 3.12.4)

# Libraries linked via full path no longer produce linker search paths.
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)

project(cosmictiger CXX C)

# Search path for CMake modules to be loaded by include() and find_package()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(HPX REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
find_package(Silo REQUIRED)
find_package(Boost REQUIRED)

enable_language(CUDA)


set(CMAKE_CXX_FLAGS "-Wno-pointer-arith -isystem /home/dmarce1/local/boost/include -DBOOST_ALL_NO_LIB -DBOOST_DISABLE_ASSERTS -O3 -std=c++17 -DNDEBUG -ffast-math -march=native")
#set(CMAKE_CXX_FLAGS "-isystem /home/dmarce1/local/boost/include -DBOOST_ALL_NO_LIB -DBOOST_DISABLE_ASSERTS -g -std=c++17 -ffast-math -march=native")

set(cosmictiger_source_files
	src/bucket.cpp
	src/check.cpp
	src/cosmos.cpp
	src/error.cpp
	src/fileio.cpp
	src/gravity.cpp
	src/gravity_cuda.cu
	src/gravity_queue.cpp
	src/hpx.cpp
	src/main.cpp
	src/options.cpp
	src/output.cpp
	src/pos_cache.cpp
	src/rand.cpp
	src/range.cpp
	src/tree.cpp
	src/tree_dir.cpp
	src/tree_client.cpp
)


set(cosmictiger_header_files
	cosmictiger/bucket.hpp	 
	cosmictiger/check.hpp	 
	cosmictiger/cosmos.hpp	 
	cosmictiger/cuda_check.hpp  
	cosmictiger/cuda_export.hpp  
	cosmictiger/defs.hpp	 
	cosmictiger/error.hpp	 
	cosmictiger/expansion.hpp	 
	cosmictiger/fileio.hpp	 
	cosmictiger/gravity.hpp
	cosmictiger/gravity_cuda.hpp
	cosmictiger/gravity_queue.hpp
 	cosmictiger/green.hpp
 	cosmictiger/hpx.hpp	
	cosmictiger/interactions.hpp	
	cosmictiger/multipole.hpp	
	cosmictiger/options.hpp	
	cosmictiger/output.hpp	
	cosmictiger/particle.hpp	
	cosmictiger/pinned_vector.hpp
	cosmictiger/pos_cache.hpp	
	cosmictiger/position.hpp	
	cosmictiger/rand.hpp	 
	cosmictiger/range.hpp	 
	cosmictiger/read_init.hpp	 
	cosmictiger/simd.hpp	 
	cosmictiger/time.hpp	 
	cosmictiger/timer.hpp	 
	cosmictiger/tree.hpp
	cosmictiger/tree_client.hpp 
	cosmictiger/tree_dir.hpp
	cosmictiger/tree_dir_impl.hpp
	cosmictiger/tree_mems.hpp
	cosmictiger/util.hpp
	cosmictiger/vect.hpp
)


add_library(hpx_lite 
	hpx_lite/src/hpx.cpp
	hpx_lite/src/id_type.cpp
	hpx_lite/src/thread.cpp
	hpx_lite/src/future.cpp
	hpx_lite/src/mutex.cpp
	hpx_lite/src/serialize.cpp
)

add_executable( cosmictiger_lite
    ${cosmictiger_source_files}
)

add_hpx_executable(
  cosmictiger
  DEPENDENCIES
    Silo::silo Boost::boost
  SOURCES
    ${cosmictiger_source_files}
  HEADERS
    ${cosmictiger_header_files}
 )



target_link_libraries (cosmictiger_lite /home/dmarce1/local/boost/lib/libboost_program_options.so tcmalloc hpx_lite siloh5 hwloc dl z hdf5)


target_include_directories(hpx_lite PUBLIC ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/hpx_lite)
target_include_directories(cosmictiger PUBLIC ${PROJECT_SOURCE_DIR})
target_compile_options(cosmictiger PUBLIC "-DUSE_HPX")
target_include_directories(cosmictiger_lite PUBLIC ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/hpx_lite)
target_compile_options(cosmictiger_lite PUBLIC "-DHPX_LITE")
target_compile_options(cosmictiger_lite PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                       --gpu-architecture=compute_75;
						-Xptxas -O3 --use_fast_math ;
						-v
                        >)
target_compile_options(cosmictiger PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                       --gpu-architecture=compute_75;
                        >)

