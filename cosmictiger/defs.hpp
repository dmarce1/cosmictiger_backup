#pragma once

#define NDIM 3
#define NCHILD 2
#define NSIBLING 6

#define USE_AVX2

#define DEFAULT_CLASS_MEMBERS(c) 		\
	c() = default; 						\
	c(const c&) = default;				\
	c(c&&) = default;					\
	c& operator=(const c&) = default;	\
	c& operator=(c&&) = default;		\
	~c() = default
