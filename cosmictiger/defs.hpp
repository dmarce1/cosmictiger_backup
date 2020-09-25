#pragma once

#define NDIM 3


#define DEFAULT_CLASS_MEMBERS(c) 		\
	c() = default; 						\
	c(const c&) = default;				\
	c(c&&) = default;					\
	c& operator=(const c&) = default;	\
	c& operator=(c&&) = default;		\
	~c() = default
