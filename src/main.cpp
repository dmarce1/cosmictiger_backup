#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>

#include <cosmictiger/rand.hpp>
#include <cosmictiger/timer.hpp>
#include <cosmictiger/tree.hpp>

#include <ctime>

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
//	srand (time(NULL));
	bucket parts;
	for (int i = 0; i < opts.problem_size; i++) {
		particle p;
		p.x = double_to_pos(rand_unit_vect());
		p.v = rand_unit_vect() * std::pow(std::max(opts.problem_size / opts.bucket_size, std::uint64_t(1)), - 1.0 / 3.0) / 8.0;
		p.out = 0;
		p.step = 0;
		p.group = 0;
		parts.insert(p);
	}
	range root_box;
	for (int dim = 0; dim < NDIM; dim++) {
		root_box.min[dim] = 0.0;
		root_box.max[dim] = 1.0;
	}
	auto root_id = hpx::new_ < tree > (hpx::find_here(), root_box, 0).get();
	auto root = tree_client(std::move(root_id), hpx::async < tree::get_ptr_action > (root_id).get());
	auto ts = timer();
	printf("Growing\n");
	root.grow(0, bucket()).get();
	printf("Grown\n");
	constexpr std::uint64_t chunk_size = 32 * 1024 * 1024;
	for (std::uint64_t i = 0; i < opts.problem_size; i += chunk_size) {
		bucket these_parts;
		for (std::uint64_t j = 0; j < chunk_size; j++) {
			if (parts.size()) {
				these_parts.insert(parts.front());
				parts.remove(parts.begin());
			} else {
				break;
			}
		}
		root.load_balance(0, 0).get();
		printf("Balanced\n");
		auto count = root.grow(0, std::move(these_parts)).get();
		printf("Counted %li parts %e s\n", count, timer() - ts);
	}
	printf("Grown\n");
	root.load_balance(0, 0).get();
	printf("Balanced\n");
	ts = timer();
	int step = 0;
	ts = timer();
	printf("Drifting\n");
	std::uint64_t cnt = root.drift(0, step++, tree_client(), root, 0.01).get();
	printf("Pruning\n");
	root.prune(0).get();
	double dtime = timer() - ts;
	double pct_drift = double(cnt) / opts.problem_size * 100;
	printf("Drift takes %e seconds %f%% drifted\n", dtime, pct_drift);
	ts = timer();
	int rc = root.verify(0).get();
	if (rc) {
		printf("%s\n", tree_verification_error(rc).c_str());
	}
	double traverse = timer() - ts;
	printf("Tree traversal takes %e seconds\n", traverse);
	auto count = root.grow(0, bucket()).get();
	printf("Counted %li parts\n", count);
	root.destroy(0).get();
	FILE *fp = fopen("data.txt", "at");
	fprintf(fp, "%li %e %e %e\n", opts.problem_size, traverse, dtime, pct_drift);
	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector < std::string > cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
	printf("exiting.2..\n");
}
#endif
