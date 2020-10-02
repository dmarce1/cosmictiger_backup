#include <cosmictiger/tree.hpp>
#include <cosmictiger/fileio.hpp>
#include <cosmictiger/hpx.hpp>
#include <cosmictiger/options.hpp>

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

void fileio_init_read(const std::string basename) {
	auto locality = hpx::get_locality_id();
	const int il = ((locality + 1) << 1) - 1;
	const int ir = ((locality + 1) << 1);
	std::vector<hpx::future<void>> futs;
	if (il < hpx_localities().size()) {
		futs.push_back(hpx::async < fileio_init_read_action > (hpx_localities()[il], basename));
	}
	if (ir < hpx_localities().size()) {
		futs.push_back(hpx::async < fileio_init_read_action > (hpx_localities()[ir], basename));
	}

	std::int32_t dummy;
	std::string filename = basename + '.' + std::to_string(locality);
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
	if (locality = 0) {
		printf("Reading %li particles\n", total_parts);
		printf("Z =             %e\n", header.redshift);
		printf("particle mass = %e\n", header.mass[1]);
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
		if (x >= 1.0 || x < 0.0) {
			printf("Particle x out of range!\n");
			abort();
		}
		if (y >= 1.0 || y < 0.0) {
			printf("Particle x out of range!\n");
			abort();
		}
		if (z >= 1.0 || z < 0.0) {
			printf("Particle x out of range!\n");
			abort();
		}
		particle part;
		part.x[0] = x;
		part.x[1] = y;
		part.x[2] = z;
		part.rung = 0;
		part.out = 0;
		part.step = 0;
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
		futs.push_back(hpx::async < fileio_insert_parts_action > (hpx_localities()[il], chunk_size));
	} else {
		futs.push_back(hpx::make_ready_future(std::uint64_t(0)));
	}
	if (ir < hpx_localities().size()) {
		futs.push_back(hpx::async < fileio_insert_parts_action > (hpx_localities()[ir], chunk_size));
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

