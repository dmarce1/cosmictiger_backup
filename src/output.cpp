#include <cosmictiger/output.hpp>
#include <cosmictiger/options.hpp>
#include <silo.h>

static mutex_type mtx;
static std::vector<output_part> parts;

void output_add_particle(const particle &p, const _4force &f) {
	std::lock_guard<mutex_type> lock(mtx);
	output_part o;
	o.x = pos_to_double(p.x);
	o.v = p.v;
	o.rung = p.rung;
	o.phi = f.phi;
	o.g = f.g;
	o.group_id = p.group;
	parts.push_back(o);
}

std::vector<output_part> gather_output();

HPX_PLAIN_ACTION (gather_output);

std::vector<output_part> gather_output() {
	std::vector < hpx::future<std::vector<output_part>> > futs;
	const auto il = ((hpx::get_locality_id() + 1) << 1) - 1;
	const auto ir = ((hpx::get_locality_id() + 1) << 1);
	if (il < hpx_localities().size()) {
		futs.push_back(hpx::async<gather_output_action>(hpx_localities()[il]));
	}
	if (ir < hpx_localities().size()) {
		futs.push_back(hpx::async<gather_output_action>(hpx_localities()[ir]));
	}
	for (int i = 0; i < futs.size(); i++) {
		const auto these_parts = futs[i].get();
		parts.insert(parts.begin(), these_parts.begin(), these_parts.end());
	}
	return std::move(parts);
}

void output_to_file(const std::string& filename, const std::vector<output_part>& parts) {

	std::thread([&] {
		static const float m = opts.particle_mass;
		std::array<std::vector<double>, NDIM> x;
		std::array<std::vector<float>, NDIM> g;
		std::array<std::vector<float>, NDIM> v;
		std::vector<float> phi;
		std::vector<int> rung;
		std::vector<std::uint64_t> group;

		printf( "SILO out\n");
		for (auto i = parts.begin(); i != parts.end(); i++) {
			for (int dim = 0; dim < NDIM; dim++) {
				x[dim].push_back(i->x[dim]);
				g[dim].push_back(i->g[dim]);
				v[dim].push_back(i->v[dim]);
			}
			if( i->group_id == DEFAULT_GROUP) {
				group.push_back(-1);
			} else {
				group.push_back(i->group_id);
			}
			rung.push_back(i->rung);
			phi.push_back(i->phi);
		}
		DBfile *db = DBCreateReal(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Meshless", DB_PDB);
		const int nparts = phi.size();
		double *coords[NDIM] = {x[0].data(), x[1].data(), x[2].data()};
		DBPutPointmesh(db, "points", NDIM, coords, nparts, DB_DOUBLE, NULL);
		for (int dim = 0; dim < NDIM; dim++) {
			std::string nm = std::string() + "v_" + char('x' + char(dim));
			DBPutPointvar1(db, nm.c_str(), "points", v[dim].data(), nparts, DB_FLOAT, NULL);
		}
		for (int dim = 0; dim < NDIM; dim++) {
			std::string nm = std::string() + "g_" + char('x' + char(dim));
			DBPutPointvar1(db, nm.c_str(), "points", g[dim].data(), nparts, DB_FLOAT, NULL);
		}
		DBPutPointvar1(db, "phi", "points", phi.data(), nparts, DB_FLOAT, NULL);
		DBPutPointvar1(db, "rung", "points", rung.data(), nparts, DB_INT, NULL);
		DBPutPointvar1(db, "group_id", "points", group.data(), nparts, DB_LONG_LONG, NULL);
		DBClose(db);}
	).join();


}
