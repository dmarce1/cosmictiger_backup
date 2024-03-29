#include <cosmictiger/tree.hpp>
#include <cosmictiger/fileio.hpp>
#include <cosmictiger/gravity.hpp>
#include <cosmictiger/hpx.hpp>
#include <cosmictiger/rand.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/cosmos.hpp>

// This header structure was copied from N-GenIC

struct io_header {
	std::uint32_t npart[6]; /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
	double mass[6]; /*!< mass[1] gives the particle mass */
	double time; /*!< time (=cosmological scale factor) of snapshot */
	double redshift; /*!< redshift of snapshot */
	std::int32_t flag_sfr; /*!< flags whether star formation is used (not available in L-Gadget2) */
	std::int32_t flag_feedback; /*!< flags whether feedback from star formation is included */
	std::uint32_t npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
	 the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
	std::int32_t flag_cooling; /*!< flags whether radiative cooling is included */
	std::int32_t num_files; /*!< determines the number of files that are used for a snapshot */
	double BoxSize; /*!< Simulation box size (in code units) */
	double Omega0; /*!< matter density */
	double OmegaLambda; /*!< vacuum energy density */
	double HubbleParam; /*!< little 'h' */
	std::int32_t flag_stellarage; /*!< flags whether the age of newly formed stars is recorded and saved */
	std::int32_t flag_metals; /*!< flags whether metal enrichment is included */
	std::int32_t hashtabsize; /*!< gives the size of the hashtable belonging to this snapshot file */
	char fill[84]; /*!< fills to 256 Bytes */
};

static bucket parts;

HPX_PLAIN_ACTION (fileio_init_read);
HPX_PLAIN_ACTION (fileio_insert_parts);

#define FREAD_ASSERT(b) { \
	if( (b) != 1) { \
		printf( "Unable to read file. source %s line %i", __FILE__, __LINE__); \
	} \
}

void fileio_init_read() {
	auto locality = hpx::get_locality_id();
	const int il = ((locality + 1) << 1) - 1;
	const int ir = ((locality + 1) << 1);
	std::vector<hpx::future<void>> futs;
	if (il < hpx_localities().size()) {
		futs.push_back(hpx::async<fileio_init_read_action>(hpx_localities()[il]));
	}
	if (ir < hpx_localities().size()) {
		futs.push_back(hpx::async<fileio_init_read_action>(hpx_localities()[ir]));
	}

	std::int32_t dummy;
	std::string filename = opts.input_file + '.' + std::to_string(locality);
	FILE *fp = fopen(filename.c_str(), "rb");
	if (!fp) {
		printf("Locality %i: Unable to load %s\n", locality, filename.c_str());
		abort();
	}
	io_header header;
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	FREAD_ASSERT(fread(&header, sizeof(header), 1, fp));
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	const std::uint64_t total_parts = std::uint64_t(header.npartTotal[1]) + (std::uint64_t(header.npartTotal[2]) << std::uint64_t(32));
	opts.problem_size = total_parts;
	opts.h = SELF_PHI * opts.soft_len * std::pow(opts.problem_size, -1.0 / 3.0);
	opts.particle_mass = header.mass[1];
	const auto Gcgs = 6.672e-8;
	const auto Hcgs = 3.2407789e-18;
	const auto mtot = opts.particle_mass * opts.problem_size;
	opts.code_to_cm = std::pow((8.0 * M_PI * Gcgs * mtot * opts.code_to_g) / (3.0 * opts.omega_m * Hcgs * Hcgs), 1.0 / 3.0);
	opts.code_to_s = opts.code_to_cm / 3e10;
	opts.H0 = Hcgs * opts.code_to_s;
	opts.G = Gcgs / pow(opts.code_to_cm, 3) * opts.code_to_g * pow(opts.code_to_s, 2);
	opts.z0 = header.redshift;
	cosmos cosmo;
	cosmo.advance_to_scalefactor(1.0/(opts.z0+1.0));
	opts.t_max = -cosmo.t;
	if (locality == 0) {
		printf("Reading %li particles\n", total_parts);
		printf("code_to_cm =    %e\n", opts.code_to_cm);
		printf("code_to_g  =    %e\n", opts.code_to_g);
		printf("code_to_s  =    %e\n", opts.code_to_s);
		printf("t_max  =        %e\n", opts.t_max);
		printf("G          =    %e\n", opts.G);
		printf("H0         =    %e\n", opts.H0);
		printf("Z =             %e\n", header.redshift);
		printf("particle mass = %e\n", header.mass[1]);
		printf("mtot =          %e\n", mtot);
		printf("hsoft =         %e\n", opts.h);
		printf("Omega_m =       %e\n", header.Omega0);
		printf("Omega_lambda =  %e\n", header.OmegaLambda);
		printf("Hubble Param =  %e\n", header.HubbleParam);
	}
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	for (int i = 0; i < header.npart[1]; i++) {
		float x, y, z;
		FREAD_ASSERT(fread(&x, sizeof(float), 1, fp));
		FREAD_ASSERT(fread(&y, sizeof(float), 1, fp));
		FREAD_ASSERT(fread(&z, sizeof(float), 1, fp));
		if (x > 1.0 || x < 0.0) {
			printf("Particle x out of range %e!\n", x);
			abort();
		}
		if (y > 1.0 || y < 0.0) {
			printf("Particle y out of range %e!\n", y);
			abort();
		}
		if (z > 1.0 || z < 0.0) {
			printf("Particle z out of range %e!\n", z);
			abort();
		}
		if (x == 1.0) {
			x = 0.0;
		}
		if (y == 1.0) {
			y = 0.0;
		}
		if (z == 1.0) {
			z = 0.0;
		}
		particle part;
		part.x[0] = x;
		part.x[1] = y;
		part.x[2] = z;
		part.rung = 0;
		part.out = 0;
//		part.out = rand1() < opts.out_pct;
		part.group = 0;
		parts.insert(part);
	}
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	const auto c0 = std::pow(1.0 / (1.0 + header.redshift), 1.5);
	for (auto &part : parts) {
		float vx, vy, vz;
		FREAD_ASSERT(fread(&vx, sizeof(float), 1, fp));
		FREAD_ASSERT(fread(&vy, sizeof(float), 1, fp));
		FREAD_ASSERT(fread(&vz, sizeof(float), 1, fp));
		part.v[0] = vx * c0;
		part.v[1] = vy * c0;
		part.v[2] = vz * c0;
	}
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	fclose(fp);
	hpx::wait_all(futs.begin(), futs.end());
}

std::uint64_t fileio_insert_parts(std::uint64_t chunk_size) {
	auto locality = hpx::get_locality_id();
	const int il = ((locality + 1) << 1) - 1;
	const int ir = ((locality + 1) << 1);
	std::vector < hpx::future < std::uint64_t >> futs;
	if (il < hpx_localities().size()) {
		futs.push_back(hpx::async<fileio_insert_parts_action>(hpx_localities()[il], chunk_size));
	} else {
		futs.push_back(hpx::make_ready_future(std::uint64_t(0)));
	}
	if (ir < hpx_localities().size()) {
		futs.push_back(hpx::async<fileio_insert_parts_action>(hpx_localities()[ir], chunk_size));
	} else {
		futs.push_back(hpx::make_ready_future(std::uint64_t(0)));
	}

	std::uint64_t read_count = 0;

	bucket these_parts;
	while (parts.size() && these_parts.size() < chunk_size) {
		these_parts.insert(parts.front());
		parts.remove(parts.begin());
		read_count++;
	}
	tree_insert_parts(std::move(these_parts));

	for (auto &fut : futs) {
		read_count += fut.get();
	}
	return read_count;
}

