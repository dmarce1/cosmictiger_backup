#include <cosmictiger/output.hpp>

mutex_type mtx;
static std::vector<output_part> parts;

void output_add_particle(const particle &p, const _4force &f) {
	std::lock_guard<mutex_type> lock(mtx);
	output_part o;
	o.x = pos_to_double(p.x);
	o.v = p.v;
	o.rung = p.rung;
	o.phi = f.phi;
	o.g = f.g;
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

void output_to_file(const std::string&) {
	auto parts = gather_output();




}
