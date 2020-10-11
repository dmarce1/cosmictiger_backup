#pragma once
#include <cosmictiger/defs.hpp>

#include <string>

class options {
public:
	bool ewald;
	bool test;
	int bucket_size;
	double code_to_cm;
	double out_pct;
	double particle_mass;
	double soft_len;
	double h;
	double sink_bias;
	double theta;
	std::uint64_t problem_size;
	std::string input_file;
	std::string config_file;

	template<class Arc>
	void serialize(Arc &arc, unsigned) {
		arc & test;
		arc & out_pct;
		arc & ewald;
		arc & bucket_size;
		arc & code_to_cm;
		arc & particle_mass;
		arc & soft_len;
		arc & h;
		arc & theta;
		arc & problem_size;
		arc & input_file;
		arc & config_file;
		arc & sink_bias;
	}

	static void set(options);
	bool process_options(int argc, char *argv[]);
};

#ifndef OPTIONS_CPP
extern options opts;
#endif
