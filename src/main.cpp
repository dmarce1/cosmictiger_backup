#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>

#include <cosmictiger/rand.hpp>
#include <cosmictiger/timer.hpp>
#include <cosmictiger/tree.hpp>
#include <cosmictiger/fileio.hpp>

#include <ctime>

int hpx_main(int argc, char *argv[]) {
	opts.process_options(argc, argv);

	printf("Creating root node\n");
	range root_box;
	for (int dim = 0; dim < NDIM; dim++) {
		root_box.min[dim] = 0.0;
		root_box.max[dim] = 1.0;
	}
	auto root_id = hpx::new_ < tree > (hpx::find_here(), 1, root_box, 0).get();
	auto root = tree_client(root_id, hpx::async < tree::get_ptr_action > (root_id).get());

	printf("Creating distributed nodes\n");
	root.grow(0, false, bucket()).get();
	printf("Balancing\n");
	root.load_balance(0, false, 0, 1).get();
	printf("Building tree directory\n");
	hpx::async < tree::build_tree_dir_action > (root_id, root).get();
	auto ts = timer();
	printf("Reading input data\n");
	fileio_init_read();
	printf("Distributing data\n");
	std::uint64_t ndistrib;
	std::uint64_t chunk_size = 1024 * 1024;
	std::uint64_t ntotal;
	do {
		int rc = root.verify(0, false).get();
		if (rc) {
			printf("%s\n", tree_verification_error(rc).c_str());
		}
		ndistrib = fileio_insert_parts(chunk_size);
		ntotal = root.grow(0, false, bucket()).get();
		printf("%li distributed\n", ntotal);
		chunk_size = std::min(2 * chunk_size, std::uint64_t(32 * 1024 * 1024));
	} while (ndistrib > 0);
	printf("%li loaded\n", ntotal);
	printf("load balancing\n");
	root.load_balance(0, false, 0, ntotal).get();
	printf("Data loaded and distributed in %e seconds\n", timer() - ts);

//	int step = 0;
//	auto dtime = timer();
//	std::uint64_t cnt = root.drift(0, false, step++, tree_client(), root, 0.01).get();
//	printf("Pruning\n");
//	root.prune(0, false).get();
//	printf("Balancing %li particles\n", opts.problem_size);
//	root.load_balance(0, false, 0, opts.problem_size).get();
//	double pct_drift = double(cnt) / opts.problem_size * 100;
//	dtime = timer() - dtime;
//	printf("Drift takes %e seconds %f%% drifted\n", dtime, pct_drift);

	root.destroy(0).get();
	return hpx::finalize();
}
//
//int hpx_main(int argc, char *argv[]) {
//	options opts;
//	opts.process_options(argc, argv);
////	bucket test;
////	for (int i = 0; i < 17; i++) {
////		particle p;
////		p.x[0] = (i+1)/ 1000.0;
////		test.insert(p);
//////		printf("%e\n", double(p.x[0]));
////	}
////	int j = 0;
////	for (auto i = test.begin(); i != test.end();) {
////		printf( "%e\n", pos_to_double(i->x[0]) );
////		if( j == 16) {
////			i = test.remove(i);
////		} else {
////			 i++;
////		}
////		j++;
////	}
////	printf( "-------\n");
////	for (auto i = test.begin(); i != test.end(); i++) {
////		printf( "%e\n", pos_to_double(i->x[0]) );
////	}
//	srand (1234);
//	bucket parts;
//	opts.problem_size = 1024;
//	for (int i = 0; i < opts.problem_size; i++) {
//		particle p;
//		p.x = double_to_pos(rand_unit_vect());
//		p.v = rand_unit_vect() * std::pow(std::max(opts.problem_size / opts.bucket_size, std::uint64_t(1)), -1.0 / 3.0) / 8.0;
////		printf( "%e %e %e\n", p.v[0],  p.v[1],  p.v[2]);
//		p.out = 0;
//		p.step = 0;
//		p.group = 0;
//		parts.insert(p);
//	}
//	range root_box;
//	for (int dim = 0; dim < NDIM; dim++) {
//		root_box.min[dim] = 0.0;
//		root_box.max[dim] = 1.0;
//	}
//	auto root_id = hpx::new_ < tree > (hpx::find_here(), 1, root_box, 0).get();
//	auto root = tree_client(std::move(root_id), hpx::async < tree::get_ptr_action > (root_id).get());
//	printf("Growing\n");
//	root.grow(0, false, bucket()).get();
//	printf("Building directory\n");
//	auto dir = hpx::async < tree::build_tree_dir_action > (root_id, root).get();
//	std::uint64_t chunk_size = 64;
//	std::uint64_t count = 1;
//	for (std::uint64_t i = 0; i < opts.problem_size; i += chunk_size) {
//		chunk_size = std::min(std::uint64_t(2) * chunk_size, std::uint64_t(32 * 1024 * 1024));
//		bucket these_parts;
//		//	printf( "%li to %li\n", i, i + chunk_size );
//		for (std::uint64_t j = 0; j < chunk_size; j++) {
//			if (parts.size()) {
//				these_parts.insert(parts.front());
//				parts.remove(parts.begin());
//			} else {
//				break;
//			}
//		}
//		root.load_balance(0, false, 0, count).get();
//		printf("Balanced\n");
//		auto ts = timer();
//		count = root.grow(0, false, std::move(these_parts)).get();
//		printf("Counted %li parts %e s\n", count, timer() - ts);
//	}
//	printf("Grown\n");
//	root.load_balance(0, false, 0, opts.problem_size).get();
//	printf("Balanced\n");
//	auto ts = timer();
//	int step = 0;
//	ts = timer();
//	printf("Drifting\n");
//	std::uint64_t cnt = root.drift(0, false, step++, tree_client(), root, 0.01).get();
//	printf("Pruning\n");
//	root.prune(0, false).get();
//	printf("Balancing\n");
//	root.load_balance(0, false, 0, opts.problem_size).get();
//	double dtime = timer() - ts;
//	double pct_drift = double(cnt) / opts.problem_size * 100;
//	printf("Drift takes %e seconds %f%% drifted\n", dtime, pct_drift);
//	ts = timer();
//	int rc = root.verify(0, false).get();
//	if (rc) {
//		printf("%s\n", tree_verification_error(rc).c_str());
//	}
//	double traverse = timer() - ts;
//	printf("Tree traversal takes %e seconds\n", traverse);
//	count = root.grow(0, false, bucket()).get();
//	printf("Counted %li parts\n", count);
//	root.destroy(0).get();
//	FILE *fp = fopen("data.txt", "at");
//	fprintf(fp, "%li %e %e %e\n", opts.problem_size, traverse, dtime, pct_drift);
//	return hpx::finalize();
//}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector < std::string > cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
	printf("exiting.2..\n");
}
#endif
