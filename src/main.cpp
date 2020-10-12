#include <cosmictiger/defs.hpp>
#include <cosmictiger/hpx.hpp>

#include <cosmictiger/cosmos.hpp>
#include <cosmictiger/error.hpp>
#include <cosmictiger/rand.hpp>
#include <cosmictiger/time.hpp>
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

void set_params(double theta, int min_rung, bool stats, double a) {
	fmm_params params;
	params.theta = theta;
	params.min_rung = min_rung;
	params.stats = stats;
	params.a = a;
	tree_set_fmm_params(params);
}

void solve_gravity(tree_client root) {
	multipole_time -= timer();
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
	root.kick_fmm(0, false, std::move(dchecklist), std::move(echecklist), std::move(L)).get();
	fmm_time += timer();
	check_cleanup();

}

int hpx_main(int argc, char *argv[]) {
	opts.process_options(argc, argv);

	set_params(opts.theta, 0, false, 1.0);
	range root_box;
	for (int dim = 0; dim < NDIM; dim++) {
		root_box.min[dim] = 0.0;
		root_box.max[dim] = 1.0;
	}
	auto root_id = hpx::new_<tree>(hpx::find_here(), 1, root_box, 0).get();
	auto root = tree_client(root_id, hpx::async<tree::get_ptr_action>(root_id).get());

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
	if (opts.test) {
		set_params(opts.theta, 0, true, 1.0);
		solve_gravity(root);
		printf("Multipoles took %e seconds\n", multipole_time);
		printf("FMM took %e seconds\n", fmm_time);
		auto output = gather_output();
		compute_error (output);
	} else {
		int i = 0;
		double t = 0.0;
		time_type itime = 0;
		constexpr int nout = 64;
		cosmos cosmo;
		cosmo.advance_to_scalefactor(1.0 / (opts.z0 + 1.0));
		int oi = 0;
		do {
			bool output = true;
			set_params(opts.theta, min_rung(itime), output, cosmo.a);
			solve_gravity(root);
			if (output) {
				auto output = gather_output();
				compute_error (output);
				output_to_file((std::string("parts.") + std::to_string(oi) + ".silo"), output);
				oi++;
			}
			auto krc = tree_kick_return();
			double dt = rung_to_dt(krc.rung);
			double a0 = cosmo.a;
			cosmo.advance_to_time(t - opts.t_max + dt);
			double a1 = cosmo.a;
			double abar = std::sqrt(1.0 / (0.5 / (a0 * a0) + 0.5 / (a1 * a1)));
			double drift_time = timer();
	//		printf( "ndrift = %i\n", root.drift(0, false, dt, abar).get().ndrift);
			root.count_children().get();
			drift_time = timer() - drift_time;
			itime = inc(itime, krc.rung);
			t = time_to_double(itime);
			printf("%i %e %e %i %i %e %e %e\n", i, dt, cosmo.a, krc.rung, min_rung(itime), multipole_time, fmm_time, drift_time);
			multipole_time = fmm_time = 0.0;
			i++;
		} while (t < opts.t_max);
	}

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
