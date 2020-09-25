#pragma once
#include <cosmictiger/defs.hpp>


#include <string>

class options {
public:
	int bucket_size;
	std::string config_file;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & bucket_size;
		arc & config_file;
	}

	static void set(options);
	bool process_options(int argc, char *argv[]);
};


#ifndef OPTIONS_CPP
extern options opts;
#endif
