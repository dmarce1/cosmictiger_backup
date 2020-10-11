#include <cosmictiger/gravity.hpp>
#include <cosmictiger/gravity_cuda.hpp>
#include <cosmictiger/green.hpp>
#include <cosmictiger/options.hpp>
#include <cosmictiger/simd.hpp>

static const auto one = simd_float(1.0);
static const auto half = simd_float(0.5);
static const simd_float eps = simd_float(std::numeric_limits<float>::min());

#include <cosmictiger/interactions.hpp>

std::uint64_t gravity_CP_direct(expansion<float> &L, const vect<position> &x, std::vector<vect<position>> y, bool do_phi) {
	if (y.size() == 0) {
		return 0;
	}
	std::uint64_t flop = 0;
	static const auto m = opts.particle_mass;
	vect<simd_int> Y, X;
	simd_float M;
	expansion<simd_float> Lacc;
	Lacc = simd_float(0);
	const auto cnt1 = y.size();
	const auto cnt2 = ((cnt1 - 1 + simd_float::size()) / simd_float::size()) * simd_float::size();
	y.resize(cnt2);
	for (int j = cnt1; j < cnt2; j++) {
		y[j] = y[cnt1 - 1];
	}

	for (int dim = 0; dim < NDIM; dim++) {
		X[dim] = int(x[dim]);
	}
	M = m;
	for (int j = 0; j < cnt1; j += simd_float::size()) {
		for (int dim = 0; dim < NDIM; dim++) {
			for (int k = 0; k < simd_float::size(); k++) {
				Y[dim][k] = y[j + k][dim];
			}
		}
		if (j + simd_float::size() > cnt1) {
			for (int k = cnt1; k < cnt2; k++) {
				M[k - j] = 0.0;
			}
		}
		vect<simd_float> dX;
		for (int dim = 0; dim < NDIM; dim++) {
			dX[dim] = simd_float(X[dim] - Y[dim]) * simd_float(POS_INV); // 3
		}
		flop += 3;
		flop += multipole_interaction(Lacc, M, dX, false, do_phi);												// 	401
	}

	for (int i = 0; i < LP; i++) {
		L[i] += Lacc[i].sum();
	}
	flop += LP * cnt1;
	y.resize(cnt1);
	return flop;
}

std::uint64_t gravity_CC_direct(expansion<float> &L, const vect<position> &x, std::vector<const multi_src*> &y, bool do_phi) {
	if (y.size() == 0) {
		return 0;
	}
	static const auto one = simd_float(1.0);
	static const auto half = simd_float(0.5);
	std::uint64_t flop = 0;
	static thread_local vect<simd_int> X, Y;
	multipole<simd_float> M;
	expansion<simd_float> Lacc;
	Lacc = simd_float(0);
	const auto cnt1 = y.size();
	const auto cnt2 = ((cnt1 - 1 + simd_float::size()) / simd_float::size()) * simd_float::size();
	y.resize(cnt2);
	for (int dim = 0; dim < NDIM; dim++) {
		X[dim] = int(x[dim]);
	}
	std::array<multi_src, simd_float::size()> ystage;
	vect<simd_float> dX;
	for (int j = 0; j < cnt1; j += simd_float::size()) {
		for (int k = 0; k < simd_float::size(); k++) {
			if (j + k < cnt1) {
				ystage[k] = *y[j + k];
			} else {
				ystage[k].m = 0.0;
				ystage[k].x = y[cnt1 - 1]->x;
			}
		}
		for (int k = 0; k < simd_float::size(); k++) {
			for (int n = 0; n < MP; n++) {
				M[n][k] = ystage[k].m[n];
			}
		}
		for (int k = 0; k < simd_float::size(); k++) {
			for (int dim = 0; dim < NDIM; dim++) {
				Y[dim][k] = ystage[k].x[dim];
			}
		}

		for (int dim = 0; dim < NDIM; dim++) {
			dX[dim] = simd_float(X[dim] - Y[dim]) * simd_float(POS_INV); // 3
		}
		flop += 3;
		flop += multipole_interaction(Lacc, M, dX, false, do_phi);												// 986
	}

	for (int i = 0; i < LP; i++) {
		L[i] += Lacc[i].sum();
	}
	flop += LP * cnt1;
	y.resize(cnt1);
	return flop;
}

std::uint64_t gravity_CC_ewald(expansion<float> &L, const vect<position> &x, std::vector<const multi_src*> &y, bool do_phi) {
	if (y.size() == 0) {
		return 0;
	}
	gravity_CC_ewald_cuda(L, x, y, do_phi);
	return 0;
}

