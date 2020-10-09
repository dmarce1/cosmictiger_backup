#define OPTIONS_CPP
#include <cosmictiger/gravity.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/hpx.hpp>
#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>

HPX_PLAIN_ACTION(options::set, set_options_action);

options opts;

void options::set(options o) {
	opts = o;
}

bool options::process_options(int argc, char *argv[]) {
	namespace po = boost::program_options;

	po::options_description command_opts("options");

	command_opts.add_options() //
	("help", "produce help message") //
	("config_file", po::value < std::string > (&config_file)->default_value(""), "configuration file") //
	("code_to_cm", po::value<double>(&code_to_cm)->default_value(1.0), "code units to centimers conversion factor") //
	("soft_len", po::value<double>(&soft_len)->default_value(0.02), "Plummer softening length in units of mean particle separation") //
	("input_file", po::value < std::string > (&input_file)->default_value(""), "base name for input files from N-GenIC") //
	("bucket_size", po::value<int>(&bucket_size)->default_value(64), "maximum number of particles on a node") //
			;

	boost::program_options::variables_map vm;
	po::store(po::parse_command_line(argc, argv, command_opts), vm);
	po::notify(vm);
	if (vm.count("help")) {
		std::cout << command_opts << "\n";
		return false;
	}
	if (!config_file.empty()) {
		std::ifstream cfg_fs { vm["config_file"].as<std::string>() };
		if (cfg_fs) {
			po::store(po::parse_config_file(cfg_fs, command_opts), vm);
		} else {
			printf("Configuration file %s not found!\n", config_file.c_str());
			return false;
		}
	}
	po::notify(vm);
	if (input_file == "") {
		input_file = std::string("../ics/4x1/ics");
	}
	const auto loc = hpx::find_all_localities();
	const auto sz = loc.size();
	std::vector<hpx::future<void>> futs;
	set(*this);
	for (int i = 1; i < sz; i++) {
		futs.push_back(hpx::async<set_options_action>(loc[i], *this));
	}
	hpx::wait_all(futs.begin(), futs.end());
#define SHOW( opt ) std::cout << std::string( #opt ) << " = " << std::to_string(opt) << '\n';
#define SHOW_STR( opt ) std::cout << std::string( #opt ) << " = " << opt << '\n';
	SHOW(bucket_size);
	SHOW_STR(config_file);
	SHOW_STR(soft_len);
//	SHOW_STR(problem_size);
	return true;
}
