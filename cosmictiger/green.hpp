#pragma once

constexpr float EWALD_REAL_N2 = 11;
constexpr float EWALD_FOUR_N2 = 7;
constexpr float EWALD_RADIUS_CUTOFF = 2.6;

#include <cosmictiger/expansion.hpp>
#include <cosmictiger/simd.hpp>

struct ewald_indices: public std::vector<vect<float>> {
	ewald_indices(int n2max, bool nozero) {
		const int nmax = sqrt(n2max) + 1;
		vect<float> h;
		for (int i = -nmax; i <= nmax; i++) {
			for (int j = -nmax; j <= nmax; j++) {
				for (int k = -nmax; k <= nmax; k++) {
					if (i * i + j * j + k * k <= n2max) {
						h[0] = i;
						h[1] = j;
						h[2] = k;
						if (!nozero || h.dot(h) > 0) {
							this->push_back(h);
						}
					}
				}
			}
		}
		std::sort(this->begin(), this->end(), [](vect<float> &a, vect<float> &b) {
			return a.dot(a) > b.dot(b);
		});
	}
};

struct periodic_parts: public std::vector<expansion<float>> {
	periodic_parts() {
		static const ewald_indices indices(EWALD_FOUR_N2, true);
		for (auto i : indices) {
			vect<float> h = i;
			const float h2 = h.dot(h);                     // 5 OP
			expansion<float> D;
			D = 0.0;
			if (h2 > 0) {
				const float c0 = 1.0 / h2 * exp(-M_PI * M_PI * h2 / 4.0);
				D() = -(1.0 / M_PI) * c0;
				for (int a = 0; a < NDIM; a++) {
					D(a) = 2.0 * h[a] * c0;
					for (int b = 0; b <= a; b++) {
						D(a, b) = 4.0 * M_PI * h[a] * h[b] * c0;
						for (int c = 0; c <= b; c++) {
							D(a, b, c) = -8.0 * M_PI * M_PI * h[a] * h[b] * h[c] * c0;
							for (int d = 0; d <= c; d++) {
								D(a, b, c, d) = -16.0 * M_PI * M_PI * M_PI * h[a] * h[b] * h[c] * h[d] * c0;
							}
						}
					}

				}
			}
			this->push_back(D);
		}
	}
};

template<class T>  // 167
CUDA_EXPORT int green_deriv_direct(expansion<T> &D, const T &d0, const T &d1, const T &d2, const T &d3, const T &d4, const vect<T> &dx) {
	T threedxadxb;
	T dxadxbdxc;
	const auto dx0dx0 = dx[0] * dx[0];
	const auto dx0dx1 = dx[0] * dx[1];
	const auto dx0dx2 = dx[0] * dx[2];
	const auto dx1dx1 = dx[1] * dx[1];
	const auto dx1dx2 = dx[1] * dx[2];
	const auto dx2dx2 = dx[2] * dx[2];
	const auto& dx1dx0 = dx0dx1;
	const auto& dx2dx0 = dx0dx2;
	const auto& dx2dx1 = dx1dx2;
	D[0] = d0;
	D[1] = dx[0] * d1;
	D[4] = dx0dx0 * d2;
	dxadxbdxc = dx0dx0 * dx[0];
	D[10] = dxadxbdxc * d3;
	D[20] = dxadxbdxc * dx[0] * d4;
	D[2] = dx[1] * d1;
	D[5] = dx1dx0 * d2;
	dxadxbdxc = dx1dx0 * dx[0];
	D[11] = dxadxbdxc * d3;
	D[21] = dxadxbdxc * dx[0] * d4;
	D[7] = dx1dx1 * d2;
	dxadxbdxc = dx1dx1 * dx[0];
	D[13] = dxadxbdxc * d3;
	D[23] = dxadxbdxc * dx[0] * d4;
	dxadxbdxc = dx1dx1 * dx[1];
	D[16] = dxadxbdxc * d3;
	D[26] = dxadxbdxc * dx[0] * d4;
	D[30] = dxadxbdxc * dx[1] * d4;
	D[3] = dx[2] * d1;
	D[6] = dx2dx0 * d2;
	dxadxbdxc = dx2dx0 * dx[0];
	D[12] = dxadxbdxc * d3;
	D[22] = dxadxbdxc * dx[0] * d4;
	D[8] = dx2dx1 * d2;
	dxadxbdxc = dx2dx1 * dx[0];
	D[14] = dxadxbdxc * d3;
	D[24] = dxadxbdxc * dx[0] * d4;
	dxadxbdxc = dx2dx1 * dx[1];
	D[17] = dxadxbdxc * d3;
	D[27] = dxadxbdxc * dx[0] * d4;
	D[31] = dxadxbdxc * dx[1] * d4;
	D[9] = dx2dx2 * d2;
	dxadxbdxc = dx2dx2 * dx[0];
	D[15] = dxadxbdxc * d3;
	D[25] = dxadxbdxc * dx[0] * d4;
	dxadxbdxc = dx2dx2 * dx[1];
	D[18] = dxadxbdxc * d3;
	D[28] = dxadxbdxc * dx[0] * d4;
	D[32] = dxadxbdxc * dx[1] * d4;
	dxadxbdxc = dx2dx2 * dx[2];
	D[19] = dxadxbdxc * d3;
	D[29] = dxadxbdxc * dx[0] * d4;
	D[33] = dxadxbdxc * dx[1] * d4;
	D[34] = dxadxbdxc * dx[2] * d4;

	const auto dx0d2 = dx[0] * d2;
	const auto dx1d2 = dx[1] * d2;
	const auto dx2d2 = dx[2] * d2;
	D[4] += d1;
	D[10] = fma(float(3), dx0d2, D[10]);
	D[20] = fma(float(6)*dx0dx0, d3, D[20]);
	D[20] = fma(float(2), d2, D[20]);
	D[20] += d2;
	D[7] += d1;
	D[16] = fma(float(3), dx1d2, D[16]);
	D[30] = fma(float(6)*dx1dx1, d3, D[30]);
	D[30] = fma(float(2), d2, D[30]);
	D[30] += d2;
	threedxadxb = float(3) * dx1dx0;
	D[13] += dx0d2;
	D[11] += dx1d2;
	D[26] = fma(threedxadxb, d3, D[26]);
	D[21] = fma(threedxadxb, d3, D[21]);
	D[23] += d2;
	D[23] = fma(dx0dx0, d3, D[23]);
	D[23] = fma(dx1dx1, d3, D[23]);
	D[9] += d1;
	D[19] = fma(float(3), dx2d2, D[19]);
	D[34] = fma(float(6)*dx2dx2, d3, D[34]);
	D[34] = fma(float(2), d2, D[34]);
	D[34] += d2;
	threedxadxb = float(3) * dx2dx0;
	D[15] += dx0d2;
	D[12] += dx2d2;
	D[29] = fma(threedxadxb, d3, D[29]);
	D[22] = fma(threedxadxb, d3, D[22]);
	D[25] += d2;
	D[25] = fma(dx0dx0, d3, D[25]);
	D[25] = fma(dx2dx2, d3, D[25]);
	threedxadxb = float(3) * dx2dx1;
	D[18] += dx1d2;
	D[17] += dx2d2;
	D[33] = fma(threedxadxb, d3, D[33]);
	D[31] = fma(threedxadxb, d3, D[31]);
	D[32] += d2;
	D[32] = fma(dx1dx1, d3, D[32]);
	D[32] = fma(dx2dx2, d3, D[32]);
	D[28] = fma(dx1dx0, d3, D[28]);
	D[24] = fma(dx2dx1, d3, D[24]);
	D[27] = fma(dx2dx0,  d3, D[27]);
	return 134;
}

template<class T>  // 576
CUDA_EXPORT int green_deriv_ewald(expansion<T> &D, const T &d0, const T &d1, const T &d2, const T &d3, const T &d4, const vect<T> &dx) {
	T threedxadxb;
	T dxadxbdxc;
	const auto dx0dx0 = dx[0] * dx[0];
	const auto dx0dx1 = dx[0] * dx[1];
	const auto dx0dx2 = dx[0] * dx[2];
	const auto dx1dx1 = dx[1] * dx[1];
	const auto dx1dx2 = dx[1] * dx[2];
	const auto dx2dx2 = dx[2] * dx[2];
	const auto& dx1dx0 = dx0dx1;
	const auto& dx2dx0 = dx0dx2;
	const auto& dx2dx1 = dx1dx2;
	D[0] += d0;
	D[1] = fma( dx[0], d1, D[1]);
	D[4] = fma( dx0dx0, d2, D[4]);
	dxadxbdxc = dx0dx0 * dx[0];
	D[10] = fma( dxadxbdxc, d3, D[10]);
	D[20] = fma( dxadxbdxc * dx[0], d4, D[20]);
	D[2] = fma( dx[1], d1, D[2]);
	D[5] = fma( dx1dx0, d2, D[5]);
	dxadxbdxc = dx1dx0 * dx[0];
	D[11] = fma( dxadxbdxc, d3, D[11]);
	D[21] = fma( dxadxbdxc * dx[0], d4, D[21]);
	D[7] = fma( dx1dx1, d2, D[7]);
	dxadxbdxc = dx1dx1 * dx[0];
	D[13] = fma( dxadxbdxc, d3, D[13]);
	D[23] = fma( dxadxbdxc * dx[0], d4, D[23]);
	dxadxbdxc = dx1dx1 * dx[1];
	D[16] = fma( dxadxbdxc, d3, D[16]);
	D[26] = fma( dxadxbdxc * dx[0], d4, D[26]);
	D[30] = fma( dxadxbdxc * dx[1], d4, D[30]);
	D[3] = fma( dx[2], d1, D[3]);
	D[6] = fma( dx2dx0, d2, D[6]);
	dxadxbdxc = dx2dx0 * dx[0];
	D[12] = fma( dxadxbdxc, d3, D[12]);
	D[22] = fma( dxadxbdxc * dx[0], d4, D[22]);
	D[8] = fma( dx2dx1, d2, D[8]);
	dxadxbdxc = dx2dx1 * dx[0];
	D[14] = fma( dxadxbdxc, d3, D[14]);
	D[24] = fma( dxadxbdxc * dx[0], d4, D[24]);
	dxadxbdxc = dx2dx1 * dx[1];
	D[17] = fma( dxadxbdxc, d3, D[17]);
	D[27] = fma( dxadxbdxc * dx[0], d4, D[27]);
	D[31] = fma( dxadxbdxc * dx[1], d4, D[31]);
	D[9] = fma( dx2dx2, d2, D[9]);
	dxadxbdxc = dx2dx2 * dx[0];
	D[15] = fma( dxadxbdxc, d3, D[15]);
	D[25] = fma( dxadxbdxc * dx[0], d4, D[25]);
	dxadxbdxc = dx2dx2 * dx[1];
	D[18] = fma( dxadxbdxc, d3, D[18]);
	D[28] = fma( dxadxbdxc * dx[0], d4, D[28]);
	D[32] = fma( dxadxbdxc * dx[1], d4, D[32]);
	dxadxbdxc = dx2dx2 * dx[2];
	D[19] = fma( dxadxbdxc, d3, D[19]);
	D[29] = fma( dxadxbdxc * dx[0], d4, D[29]);
	D[33] = fma( dxadxbdxc * dx[1], d4, D[33]);
	D[34] = fma( dxadxbdxc * dx[2], d4, D[34]);

	const auto dx0d2 = dx[0] * d2;
	const auto dx1d2 = dx[1] * d2;
	const auto dx2d2 = dx[2] * d2;
	D[4] += d1;
	D[10] = fma(float(3), dx0d2, D[10]);
	D[20] = fma(float(6)*dx0dx0, d3, D[20]);
	D[20] = fma(float(2), d2, D[20]);
	D[20] += d2;
	D[7] += d1;
	D[16] = fma(float(3), dx1d2, D[16]);
	D[30] = fma(float(6)*dx1dx1, d3, D[30]);
	D[30] = fma(float(2), d2, D[30]);
	D[30] += d2;
	threedxadxb = float(3) * dx1dx0;
	D[13] += dx0d2;
	D[11] += dx1d2;
	D[26] = fma(threedxadxb, d3, D[26]);
	D[21] = fma(threedxadxb, d3, D[21]);
	D[23] += d2;
	D[23] = fma(dx0dx0, d3, D[23]);
	D[23] = fma(dx1dx1, d3, D[23]);
	D[9] += d1;
	D[19] = fma(float(3), dx2d2, D[19]);
	D[34] = fma(float(6)*dx2dx2, d3, D[34]);
	D[34] = fma(float(2), d2, D[34]);
	D[34] += d2;
	threedxadxb = float(3) * dx2dx0;
	D[15] += dx0d2;
	D[12] += dx2d2;
	D[29] = fma(threedxadxb, d3, D[29]);
	D[22] = fma(threedxadxb, d3, D[22]);
	D[25] += d2;
	D[25] = fma(dx0dx0, d3, D[25]);
	D[25] = fma(dx2dx2, d3, D[25]);
	threedxadxb = float(3) * dx2dx1;
	D[18] += dx1d2;
	D[17] += dx2d2;
	D[33] = fma(threedxadxb, d3, D[33]);
	D[31] = fma(threedxadxb, d3, D[31]);
	D[32] += d2;
	D[32] = fma(dx1dx1, d3, D[32]);
	D[32] = fma(dx2dx2, d3, D[32]);
	D[28] = fma(dx1dx0, d3, D[28]);
	D[24] = fma(dx2dx1, d3, D[24]);
	D[27] = fma(dx2dx0,  d3, D[27]);
	return 134 + 35;
}

template<class T>
CUDA_EXPORT inline int green_direct(expansion<T> &D, const vect<T> &dX) {
	static const T r0 = 1.0e-9;
//	static const T H = options::get().soft_len;
	static const T nthree(-3.0);
	static const T nfive(-5.0);
	static const T nseven(-7.0);
	const T r2 = dX.dot(dX);				// 5
	const T r = sqrt(r2);					// 7
	const T rinv = (r > r0) / max(r, r0);	// 3
	const T r2inv = rinv * rinv;			// 1
	const T d0 = -rinv;						// 1
	const T d1 = -d0 * r2inv;				// 2
	const T d2 = nthree * d1 * r2inv;		// 2
	const T d3 = nfive * d2 * r2inv;		// 2
	const T d4 = nseven * d3 * r2inv;		// 2
	return 25 + green_deriv_direct(D, d0, d1, d2, d3, d4, dX);
}
#ifdef __CUDA_ARCH__

__device__ const cuda_ewald_const& cuda_get_const();

CUDA_EXPORT int green_ewald(expansion<float>& D, const vect<float> &X) {
	int flop = 0;
	const auto &cuda_const = cuda_get_const();
	const auto &four_indices = cuda_const.four_indices;
	const auto &real_indices = cuda_const.real_indices;
	const auto &hparts = cuda_const.periodic_parts;
	static const float three(3.0);
	const float fouroversqrtpi(4.0 / sqrt(M_PI));
	static const float one(1.0);
	static const float two(2.0);
	static const float nthree(-3.0);
	static const float nfour(-4.0);
	static const float nfive(-5.0);
	static const float nseven(-7.0);
	static const float neight(-8.0);
	static const float p(0.3275911);
	static const float a1(0.254829592);
	static const float a2(-0.284496736);
	static const float a3(1.421413741);
	static const float a4(-1.453152027);
	static const float a5(1.061405429);
	static const float rcut(1.0e-6);
	const float r = abs(X);				// 5
	const float zmask = r > rcut;		// 1
	expansion<float>& Dreal = D;
	expansion<float> Dfour;
	Dreal = 0.0;
	Dfour = 0.0;
	flop += 6;
	for (auto n : real_indices) {
		const vect<float> dx = X - vect<float>(n);				// 3
		const float r2 = dx.dot(dx);							// 5
		flop += 8;
		if (r2 < (EWALD_RADIUS_CUTOFF * EWALD_RADIUS_CUTOFF)) {	// 1
			const float r = sqrt(r2);					// 1
			const float cmask = one - (n.dot(n) > 0.0); // 7
			const float mask = (one - (one - zmask) * cmask); // 3
			const float rinv = mask / max(r, rcut);		// 2
			const float r2inv = rinv * rinv;			// 1
			const float r3inv = r2inv * rinv;			// 1
			const float t1 = float(1) / (float(1) + p * two * r); 	//4
			const float t2 = t1 * t1;								// 1
			const float t3 = t2 * t1;								// 1
			const float t4 = t2 * t2;								// 1
			const float t5 = t2 * t3;								// 1
			const float exp0 = expf(nfour * r2);					// 26
			const float erfc0 = (a1 * t1 + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5) * exp0; 			// 10
			const float expfactor = fouroversqrtpi * r * exp0; 	// 2
			const float e1 = expfactor * r3inv;						// 1
			const float e2 = neight * e1;							// 1
			const float e3 = neight * e2;							// 1
			const float e4 = neight * e3;							// 1
			const float d0 = -erfc0 * rinv;							// 2
			const float d1 = fma(-d0,  r2inv, e1);					// 3
			const float d2 = fma(nthree*d1, r2inv,  e2);			// 3
			const float d3 = fma(nfive*d2, r2inv,  e3);				// 3
			const float d4 = fma(nseven*d3, r2inv,  e4);			// 3
			flop += 74 + green_deriv_ewald(Dreal, d0, d1, d2, d3, d4, dx);
		}
	}
	static const float twopi = 2.0 * M_PI;

	for (int i = 0; i < EWALD_NFOUR; i++) {
		const auto &h = four_indices[i];
		const auto &hpart = hparts[i];
//		printf( "H = %e %e %e\n", h[0], h[1], h[2]);
		const float h2 = h.dot(h);
		const float hdotx = h.dot(X);
		float co;
		float so;
		sincosf(twopi * hdotx, &so, &co);				// 35
		Dfour[0] = fma(hpart[0], co, Dfour[0]);			// 2 * 35
		Dfour[1] = fma(hpart[1], so, Dfour[1]);
		Dfour[2] = fma(hpart[2], so, Dfour[2]);
		Dfour[3] = fma(hpart[3], so, Dfour[3]);
		Dfour[4] = fma(hpart[4], co, Dfour[4]);
		Dfour[5] = fma(hpart[5], co, Dfour[5]);
		Dfour[6] = fma(hpart[6], co, Dfour[6]);
		Dfour[7] = fma(hpart[7], co, Dfour[7]);
		Dfour[8] = fma(hpart[8], co, Dfour[8]);
		Dfour[9] = fma(hpart[9], co, Dfour[9]);
		Dfour[10] = fma(hpart[10], so, Dfour[10]);
		Dfour[11] = fma(hpart[11], so, Dfour[11]);
		Dfour[12] = fma(hpart[12], so, Dfour[12]);
		Dfour[13] = fma(hpart[13], so, Dfour[13]);
		Dfour[14] = fma(hpart[14], so, Dfour[14]);
		Dfour[15] = fma(hpart[15], so, Dfour[15]);
		Dfour[16] = fma(hpart[16], so, Dfour[16]);
		Dfour[17] = fma(hpart[17], so, Dfour[17]);
		Dfour[18] = fma(hpart[18], so, Dfour[18]);
		Dfour[19] = fma(hpart[19], so, Dfour[19]);
		Dfour[20] = fma(hpart[20], co, Dfour[20]);
		Dfour[21] = fma(hpart[21], co, Dfour[21]);
		Dfour[22] = fma(hpart[22], co, Dfour[22]);
		Dfour[23] = fma(hpart[23], co, Dfour[23]);
		Dfour[24] = fma(hpart[24], co, Dfour[24]);
		Dfour[25] = fma(hpart[25], co, Dfour[25]);
		Dfour[26] = fma(hpart[26], co, Dfour[26]);
		Dfour[27] = fma(hpart[27], co, Dfour[27]);
		Dfour[28] = fma(hpart[28], co, Dfour[28]);
		Dfour[30] = fma(hpart[30], co, Dfour[30]);
		Dfour[29] = fma(hpart[29], co, Dfour[29]);
		Dfour[31] = fma(hpart[31], co, Dfour[31]);
		Dfour[32] = fma(hpart[32], co, Dfour[32]);
		Dfour[33] = fma(hpart[33], co, Dfour[33]);
		Dfour[34] = fma(hpart[34], co, Dfour[34]);
		flop += 105;
	}
	for (int i = 0; i < LP; i++) {
		Dreal[i] += Dfour[i];
	}
	flop += LP;

	expansion<float> D1;
	flop += green_direct(D1, X);
	D() = (M_PI / 4.0) + D();
	for( int i = 0; i < LP; i++) {
		D[i] = fma(-zmask, D1[i], D[i]);
	}
	flop += 72;
	return flop;
}

#else
#ifndef __CUDACC__

template<class T>
inline int green_ewald(expansion<T> &rcD, const vect<T> &X) {		// 251176
	static const periodic_parts periodic;
	expansion<T> D;
	D = 0.0;
	vect<T> n;
	vect<float> h;
	static const ewald_indices indices_real(EWALD_REAL_N2, false);
	static const ewald_indices indices_four(EWALD_FOUR_N2, true);
//	printf( "%i %i\n", indices_real.size(), indices_four.size());
	static const T three(3.0);
	static const T fouroversqrtpi(4.0 / sqrt(M_PI));
	static const T two(2.0);
	static const T eight(8.0);
	static const T fifteen(15.0);
	static const T thirtyfive(35.0);
	static const T fourty(40.0);
	static const T fiftysix(56.0);
	static const T sixtyfour(64.0);
	static const T onehundredfive(105.0);
	static const float rcut(1.0e-6);

//	printf("%i %i\n", indices_real.size(), indices_four.size());
	const T r = abs(X);															// 5
	const simd_float zmask = r > rcut;											// 2
	for (int i = 0; i < indices_real.size(); i++) {			// 739 * 305 		// 225395
		h = indices_real[i];
		n = h;
		const vect<T> dx = X - n;				// 3
		const T r2 = dx.dot(dx);				// 5
		const T r = sqrt(r2);					// 7
		const T cmask = T(1) - (n.dot(n) > 0.0);
		const T mask = (T(1) - (T(1) - zmask) * cmask) * (r < EWALD_RADIUS_CUTOFF);
		const T rinv = mask / max(r, rcut);		// 36
		const T r2inv = rinv * rinv;			// 1
		const T r3inv = r2inv * rinv;			// 1
		T expfac;
		const T erfc = erfcexp(two * r, &expfac);			// 76
		const T expfactor = fouroversqrtpi * r * expfac; 	// 2
		const T e1 = expfactor * r3inv;
		const T e2 = -T(8) * e1;
		const T e3 = -T(8) * e2;
		const T e4 = -T(8) * e3;
		const T d0 = -erfc * rinv;							// 2
		const T d1 = fma(-d0, r2inv, e1);			// 2
		const T d2 = fma(-T(3)*d1, r2inv, e2);			// 2
		const T d3 = fma(-T(5)*d2, r2inv, e3);			// 2
		const T d4 = fma(-T(7)*d3, r2inv, e4);			// 2
		green_deriv_ewald(D, d0, d1, d2, d3, d4, dx);			// 576
	}
	for (int i = 0; i < indices_four.size(); i++) {		// 207 * 123 = 			// 25461
		h = indices_four[i];
		const auto H = periodic[i];
		T hdotdx = X[0] * h[0];		// 1
		for (int a = 1; a < NDIM; a++) {
			hdotdx += X[a] * h[a];								// 4
		}
		static const T twopi = 2.0 * M_PI;
		const T omega = twopi * hdotdx;							// 1
		T co, si;
		sincos(omega, &si, &co);								// 25
		D() += simd_float(H()) * co;										// 5
		for (int a = 0; a < NDIM; a++) {
			D(a) += simd_float(H(a)) * si;									// 15
			for (int b = 0; b <= a; b++) {
				D(a, b) += simd_float(H(a, b)) * co;						// 30
				for (int c = 0; c <= b; c++) {
					D(a, b, c) += simd_float(H(a, b, c)) * si;				// 50
					for (int d = 0; d <= c; d++) {
						D(a, b, c, d) += simd_float(H(a, b, c, d)) * co; 	// 75
					}
				}

			}
		}
	}
	for (int i = 0; i < LP; i++) {
		rcD[i] = D[i];																	// 70
	}
	expansion<T> D1;
	green_direct(D1, X);													// 167
	rcD() = T(M_PI / 4.0) + rcD();
	for (int i = 0; i < LP; i++) {
		rcD[i] = fma(-zmask, D1[i], rcD[i]);
	}
	return 251176;
}
#endif
#endif
