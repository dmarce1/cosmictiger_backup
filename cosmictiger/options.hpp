#pragma once
#include <cosmictiger/defs.hpp>


#include <string>

class options {
public:
	int bucket_size;
	double code_to_cm;
	double particle_mass;
	double soft_len;
	double h;
	std::uint64_t problem_size;
	std::string input_file;
	std::string config_file;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & h;
		arc & code_to_cm;
		arc & bucket_size;
		arc & input_file;
		arc & particle_mass;
		arc & config_file;
	}

	static void set(options);
	bool process_options(int argc, char *argv[]);
};


#ifndef OPTIONS_CPP
extern options opts;
#endif
