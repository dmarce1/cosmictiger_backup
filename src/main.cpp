#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>

#include <cosmictiger/rand.hpp>
#include <cosmictiger/tree.hpp>

int hpx_main(int argc, char *argv[]) {
	options opts;
	opts.process_options(argc, argv);
	bucket parts;
	for (int i = 0; i < opts.problem_size; i++) {
		particle p;
		p.x = double_to_pos(rand_unit_vect());
		p.v = rand_unit_vect() * 1.0e-3;
		p.out = 0;
		p.step = 0;
		p.group = 0;
		p.dt = 0.0;
		parts.insert(p);
	}
	printf("Bucket size = %i\n", parts.size());
	tree_client root = hpx::new_ < tree > (hpx::find_here(), 1).get();
	printf("Growing\n");
	auto count = root.grow(0, std::move(parts)).get();
	printf( "Counted %li parts\n", count);
	printf("Grown\n");
	int rc = root.verify(0).get();
	if (rc) {
		printf("%s\n", tree_verification_error(rc).c_str());
	}

	int step = 0;
	root.drift(0, step++, 0.01).get();
	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector < std::string > cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
}
#endif
