#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>


int hpx_main(int argc, char *argv[]) {
	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
}

#endif
