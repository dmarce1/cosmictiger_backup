#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>

#include <cosmictiger/error.hpp>
#include <cosmictiger/rand.hpp>
#include <cosmictiger/timer.hpp>
#include <cosmictiger/tree.hpp>
#include <cosmictiger/fileio.hpp>
#include <cosmictiger/output.hpp>

#include <ctime>

double multipole_time = 0.0;
double fmm_time = 0.0;

void yield_to_hpx() {
	hpx::this_thread::yield();
}

void set_params(double theta, int min_rung, bool stats) {
	fmm_params params;
	params.theta = theta;
	params.min_rung = min_rung;
	params.stats = stats;
	tree_set_fmm_params(params);
}

void solve_gravity(tree_client root, double theta, int min_rung, bool stats) {
	set_params(theta, min_rung, stats);
	multipole_time -= timer();
	printf("multipoles\n");
	root.compute_multipoles(0, false, -1, bucket(), 0).get();
	multipole_time += timer();
	fmm_time -= timer();
	check_item root_check;
	check_info root_info;
	root_info = root.get_check_info();
	root_check.opened = false;
	root_check.info = &root_info;
	std::vector<check_item> echecklist(1, root_check);
	std::vector<check_item> dchecklist(1, root_check);
	expansion_src L;
	L.l = 0.0;
	L.x[0] = L.x[1] = L.x[2] = 0.5;
	printf("fmm\n");
	root.kick_fmm(0, false, std::move(dchecklist), std::move(echecklist), std::move(L)).get();
	fmm_time += timer();
	printf("Cleaning up\n");
	check_cleanup();

}

int hpx_main(int argc, char *argv[]) {
	opts.process_options(argc, argv);
	set_params(opts.theta, 0, false);
	range root_box;
	for (int dim = 0; dim < NDIM; dim++) {
		root_box.min[dim] = 0.0;
		root_box.max[dim] = 1.0;
	}
	auto root_id = hpx::new_<tree>(hpx::find_here(), 1, root_box, 0).get();
	auto root = tree_client(root_id, hpx::async<tree::get_ptr_action>(root_id).get());

	if (opts.test) {
		root.grow(0, false, bucket(), true).get();
		auto tstat = root.load_balance(0, false, 0, 1).get();
		root.grow(0, false, bucket()).get();
		hpx::async<tree::build_tree_dir_action>(root_id, root).get();
		fileio_init_read();
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
			chunk_size = std::min(2 * chunk_size, std::uint64_t(32 * 1024 * 1024));
		} while (ndistrib > 0);
		printf("%li loaded\n", ntotal);
		tstat = root.load_balance(0, false, 0, ntotal).get();

		solve_gravity(root, opts.theta, 0, true);
		printf("Multipoles took %e seconds\n", multipole_time);
		printf("FMM took %e seconds\n", fmm_time);
		auto output = gather_output();
		compute_error(output);
		output_to_file("parts.silo", output);
	}

	//
//	int rc = root.verify(0, false).get();
//	if (rc) {
//		printf("%s\n", tree_verification_error(rc).c_str());
//	}

//
////	double btime = 0.0;
//	for (int step = 0; step < 10; step++) {
//		auto dtime = timer();
//		auto dr = root.drift(0, false, 0.009).get();
//		printf( "D: %e\n", timer() - dtime);
//		dtime = timer();
////		int cnt = 0;
//		auto cnt = root.count_children().get();
//		printf( "C: %e\n", timer() - dtime);
//		dtime = timer();
//		root.compute_multipoles(0, false, -1, std::move(dr.parts), 0);
//		printf( "M: %e\n", timer() - dtime);
//		dtime = timer();
////		int rc = root.verify(0, false).get();
////		if (rc) {
////			printf("%s\n", tree_verification_error(rc).c_str());
////		}
//		printf("Drift takes %e seconds %i\n", timer() - dtime, cnt);
//	}
////	printf("%e\n", btime);
	printf("Destroying tree\n");
	root.destroy(0).get();
	printf("exiting\n");
	return hpx::finalize();
}

#ifndef HPX_LITE
int main(int argc, char *argv[]) {

	std::vector < std::string > cfg = { "hpx.commandline.allow_unknown=1" };

	hpx::init(argc, argv, cfg);
	printf("exiting.2..\n");
}
#endif
