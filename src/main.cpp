#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>

#include <cosmictiger/rand.hpp>
#include <cosmictiger/timer.hpp>
#include <cosmictiger/tree.hpp>

int hpx_main(int argc, char *argv[]) {
	options opts;
	opts.process_options(argc, argv);
//	bucket test;
//	for (int i = 0; i < 17; i++) {
//		particle p;
//		p.x[0] = (i+1)/ 1000.0;
//		test.insert(p);
////		printf("%e\n", double(p.x[0]));
//	}
//	int j = 0;
//	for (auto i = test.begin(); i != test.end();) {
//		printf( "%e\n", pos_to_double(i->x[0]) );
//		if( j == 16) {
//			i = test.remove(i);
//		} else {
//			 i++;
//		}
//		j++;
//	}
//	printf( "-------\n");
//	for (auto i = test.begin(); i != test.end(); i++) {
//		printf( "%e\n", pos_to_double(i->x[0]) );
//	}

	bucket parts;
	for (int i = 0; i < opts.problem_size; i++) {
		particle p;
		p.x = double_to_pos(rand_unit_vect());
		p.v = rand_unit_vect() * 1.0e-3;
		p.out = 0;
		p.step = 0;
		p.group = 0;
		parts.insert(p);
	}
	printf("Bucket size = %i\n", parts.size());
	auto root_id = hpx::new_ < tree > (hpx::find_here(), 1).get();
	auto root = tree_client(std::move(root_id), hpx::async < tree::get_ptr_action > (root_id).get());
	auto ts = timer();
	printf("Growing\n");
	root.grow(0, bucket()).get();
	printf("Grown\n");
	root.load_balance(0, 0).get();
	printf( "Balanced\n");
	auto count = root.grow(0, std::move(parts)).get();
	printf("Counted %li parts %e s\n", count, timer() - ts);
	printf("Grown\n");
	root.load_balance(0, 0).get();
	printf( "Balanced\n");
	ts = timer();
	int step = 0;
	ts = timer();
	root.drift(0, step++, tree_client(), root, 0.01).get();
	root.prune(0).get();
	printf("Drift takes %e seconds\n", timer() - ts);
	int rc = root.verify(0).get();
	if (rc) {
		printf("%s\n", tree_verification_error(rc).c_str());
	}
	printf("Tree traversal takes %e seconds\n", timer() - ts);
	root.destroy(0).get();
	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector < std::string > cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
	printf("exiting.2..\n");
}
#endif
