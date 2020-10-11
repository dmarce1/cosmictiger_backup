#include <cosmictiger/options.hpp>
#include <cosmictiger/cuda_export.hpp>
#include <cosmictiger/cuda_check.hpp>
#include <cosmictiger/gravity_cuda.hpp>
#include <cosmictiger/green.hpp>
#include <cosmictiger/interactions.hpp>

#include <stack>
#include <atomic>
void yield_to_hpx();

__device__ __constant__ cuda_ewald_const cuda_ewald;

__device__ const cuda_ewald_const& cuda_get_const() {
	return cuda_ewald;
}

double *flop_ptr;

double cuda_reset_flop() {
	double result;
	double zero = 0.0;
	CUDA_CHECK(cudaMemcpy(&result, flop_ptr, sizeof(double), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(flop_ptr, &zero, sizeof(double), cudaMemcpyHostToDevice));
	return result;
}

void cuda_init() {
	static std::atomic<int> lock(0);
	static bool init = false;
	while (lock++ != 0) {
		lock--;
	}
	if (!init) {
		static const float efs[LP + 1] = { 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 5.00000000e-01, 1.00000000e+00, 1.00000000e+00,
				5.00000000e-01, 1.00000000e+00, 5.00000000e-01, 1.66666667e-01, 5.00000000e-01, 5.00000000e-01, 5.00000000e-01, 1.00000000e+00, 5.00000000e-01,
				1.66666667e-01, 5.00000000e-01, 5.00000000e-01, 1.66666667e-01, 4.16666667e-02, 1.66666667e-01, 1.66666667e-01, 2.50000000e-01, 5.00000000e-01,
				2.50000000e-01, 1.66666667e-01, 5.00000000e-01, 5.00000000e-01, 1.66666667e-01, 4.16666667e-02, 1.66666667e-01, 2.50000000e-01, 1.66666667e-01,
				4.16666667e-02, 0.0 };

		static cuda_ewald_const c;
		const ewald_indices indices_real(EWALD_REAL_N2, false);
		const ewald_indices indices_four(EWALD_FOUR_N2, true);
		const periodic_parts periodic;
		for (int i = 0; i < indices_real.size(); i++) {
			c.real_indices[i] = indices_real[i];
		}
		for (int i = 0; i < indices_four.size(); i++) {
			c.four_indices[i] = indices_four[i];
			c.periodic_parts[i] = periodic[i];
		}
		for (int i = 0; i < LP; i++) {
			c.exp_factors[i] = efs[i];
		}

		CUDA_CHECK(cudaMemcpyToSymbol(cuda_ewald, &c, sizeof(cuda_ewald_const)));
		init = true;
		CUDA_CHECK(cudaMallocHost((void** )&flop_ptr, sizeof(double)));
		cuda_reset_flop();
//		CUDA_CHECK(cudaThreadSetLimit(cudaLimitStackSize, 2048));
	}
	lock--;
}

#define WARPSIZE 32
#define CCSIZE 32

#define WORKSIZE 256
#define PCWORKSIZE 96
#define NODESIZE 64
#define NWARP (WORKSIZE/WARPSIZE)
#define PCNWARP (PCWORKSIZE/WARPSIZE)
#define WARPSIZE 32

#include <cstdint>

#define TILESIZE 512

__global__ void gravity_ewald_direct_kernel(_4force *f, vect<double> *x, vect<double> *y, int xsize, int ysize, float m, float h) {
	const int i = threadIdx.x;
	const auto &cuda_const = cuda_get_const();
	const auto &four_indices = cuda_const.four_indices;
	const auto &real_indices = cuda_const.real_indices;
	const auto &hparts = cuda_const.periodic_parts;
	const float fouroversqrtpi(4.0 / sqrt(M_PI));
	static const float one(1.0);
	static const float two(2.0);
	static const float nfour(-4.0);
	static const float a1(0.254829592);
	static const float a2(-0.284496736);
	static const float a3(1.421413741);
	static const float a4(-1.453152027);
	static const float a5(1.061405429);
	static const float rcut(1.0e-6);
	static const float twopi = 2.0 * M_PI;
	static const float p(0.3275911);
	const float hinv = 1.0 / h;
	const float h3inv = hinv * hinv * hinv;

	for (int j = 0; j < xsize; j += TILESIZE) {
		const int k = j + i;
		if (k < xsize) {
			f[k].phi = m * SELF_PHI / h;
			f[k].g = vect<float>(0.0);
			for (int l = 0; l < ysize; l++) {
				vect<float> X;
				for (int dim = 0; dim < NDIM; dim++) {
					const auto dx = x[k][dim] - y[l][dim];
					X[dim] = float(copysign(min(abs(dx),double(1.0) - abs(dx)), dx * (double(0.5) - abs(dx))));
				}
				const float r = abs(X);
				if (r > h) {
					_4force freal;
					_4force ffour;
					freal.g = ffour.g = vect<float>(0.0);
					freal.phi = ffour.phi = 0.0;
//					for (auto n : real_indices) {
//						const vect<float> dx = X - vect<float>(n);				// 3
//						const float r2 = dx.dot(dx);							// 5
//						if (r2 < (EWALD_RADIUS_CUTOFF * EWALD_RADIUS_CUTOFF)) {	// 1
//							const float r = sqrt(r2);					// 1
//							const float cmask = one - (n.dot(n) > 0.0); // 7
//							const float rinv = one / r;		// 2
//							const float r2inv = rinv * rinv;			// 1
//							const float r3inv = r2inv * rinv;			// 1
//							const float t1 = float(1) / (float(1) + p * two * r); 	//4
//							const float t2 = t1 * t1;								// 1
//							const float t3 = t2 * t1;								// 1
//							const float t4 = t2 * t2;								// 1
//							const float t5 = t2 * t3;								// 1
//							const float exp0 = expf(nfour * r2);					// 26
//							const float erfc0 = (a1 * t1 + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5) * exp0; 			// 10
//							const float expfactor = fouroversqrtpi * r * exp0; 	// 2
//							const float e1 = expfactor * r3inv;						// 1
//							const float d0 = -erfc0 * rinv;							// 2
//							const float d1 = fma(-d0, r2inv, e1);					// 3
//							freal.phi += d0;
//							freal.g -= dx * d1;
//						}
//					}
//					for (int n = 0; n < EWALD_NFOUR; n++) {
//						const auto &h = four_indices[n];
//						const auto &hpart = hparts[n];
//						const float h2 = h.dot(h);
//						const float hdotx = h.dot(X);
//						float co;
//						float so;
//						sincosf(twopi * hdotx, &so, &co);
//						ffour.phi += hpart() * co;
//						for (int dim = 0; dim < NDIM; dim++) {
//							ffour.g[dim] -= hpart(dim) * so;
//						}
//					}
					f[k].phi -= m / r;
					f[k].g -= X * m / (r * r * r);
//					f[k].phi += (ffour.phi + freal.phi) * m;
					//					f[k].g += (ffour.g + freal.g) * m;
				} else {
					const float rinv = 1.0 / r;
					const float rinv3 = rinv * rinv;
					float p, f0;
					if (r > 0.5 * h) {
						const float roh = min(r * hinv, 1.0);                         // 2
						const float roh2 = roh * roh;                         // 1
						const float roh3 = roh2 * roh;                         // 1
						f0 = float(-32.0 / 3.0);
						f0 = fma(f0, roh, float(+192.0 / 5.0));                         // 2
						f0 = fma(f0, roh, float(-48.0));                         // 2
						f0 = fma(f0, roh, float(+64.0 / 3.0));                         // 2
						f0 = fma(f0, roh3, float(-1.0 / 15.0));                         // 2
						f0 *= rinv3;                         // 1
						p = float(+32.0 / 15.0);
						p = fma(p, roh, float(-48.0 / 5.0));                                 // 2
						p = fma(p, roh, float(+16.0));                                 // 2
						p = fma(p, roh, float(-32.0 / 3.0));                                 // 2
						p = fma(p, roh2, float(+16.0 / 5.0));                                 // 2
						p = fma(p, roh, float(-1.0 / 15.0));                                 // 2
						p *= rinv;                                 // 1
					} else {
						const float roh = min(r * hinv, 1.0);                           // 2
						const float roh2 = roh * roh;                           // 1
						f0 = float(+32.0);
						f0 = fma(f0, roh, float(-192.0 / 5.0));                           // 2
						f0 = fma(f0, roh2, float(+32.0 / 3.0));                           // 2
						f0 *= h3inv;                           // 1
						p = float(-32.0 / 5.0);
						p = fma(p, roh, float(+48.0 / 5.0));							// 2
						p = fma(p, roh2, float(-16.0 / 3.0));							// 2
						p = fma(p, roh2, float(+14.0 / 5.0));							// 2
						p *= hinv;							// 1
					}
					const auto dXM = X * m;								// 3
					f[k].g -= dXM * f0;
					f[k].phi -= p * m;
				}
			}
		}
	}
}

void gravity_ewald_direct(std::vector<_4force> &f, const std::vector<vect<double>> x, const std::vector<vect<double>> &y) {
	_4force *fdev;
	vect<double> *xdev;
	vect<double> *ydev;
	const auto xbytes = x.size() * sizeof(vect<double> );
	const auto ybytes = y.size() * sizeof(vect<double> );
	const auto fbytes = f.size() * sizeof(_4force);
	CUDA_CHECK(cudaMalloc((void** ) &fdev, fbytes));
	CUDA_CHECK(cudaMalloc((void** ) &xdev, xbytes));
	CUDA_CHECK(cudaMalloc((void** ) &ydev, ybytes));
	CUDA_CHECK(cudaMemcpy(xdev, x.data(), xbytes, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(ydev, y.data(), ybytes, cudaMemcpyHostToDevice));
	/**/gravity_ewald_direct_kernel<<<1,TILESIZE>>>(fdev, xdev, ydev, x.size(), y.size(), opts.particle_mass, opts.h);

	CUDA_CHECK(cudaMemcpy(f.data(), fdev, fbytes, cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaFree(fdev));
	CUDA_CHECK(cudaFree(xdev));
	CUDA_CHECK(cudaFree(ydev));
}

__global__ void CC_ewald_kernel(expansion<float> *lptr, const vect<position> X, const multi_src *y, int ysize, bool do_phi, double *flop_ptr) {

	int l = threadIdx.x + blockDim.x * blockIdx.x;
	int n = threadIdx.x;
	int tb_size = blockDim.x;
	auto &L = *lptr;

	__shared__ expansion<float>
	Lacc[CCSIZE];
	__shared__ std::uint64_t
	flop[CCSIZE];
	flop[n] = 0;
	for (int i = 0; i < LP; i++) {
		Lacc[n][i] = 0.0;
	}
	for (int yi = l; yi < ysize; yi += tb_size * gridDim.x) {
		if (yi < ysize) {
			vect<float> dX;
			for (int dim = 0; dim < NDIM; dim++) {
				dX[dim] = float(int(X[dim]) - int(y[yi].x[dim])) * float(POS_INV); 		// 3
			}
			flop[n] += 3 + multipole_interaction(Lacc[n], y[yi].m, dX, true, do_phi);
		}
	}
	for (int N = tb_size / 2; N > 0; N >>= 1) {
		if (n < N) {
			for (int i = 0; i < LP; i++) {
				Lacc[n][i] += Lacc[n + N][i];
			}
			flop[n] += LP;
		}
	}
	if (n == 0) {
		for (int i = 0; i < LP; i++) {
			atomicAdd(&L[i], Lacc[0][i]);
		}
		flop[n] += LP;
	}
	for (int N = tb_size / 2; N > 0; N >>= 1) {
		if (n < N) {
			flop[n] += flop[n + N];
		}
	}
	if (n == 0) {
		atomicAdd(flop_ptr, flop[0]);
	}
}

struct cuda_context_ewald {
	int ysize;
	cudaStream_t stream;
	expansion<float> *L;
	multi_src *y;
	expansion<float> *Lp;
	multi_src *yp;
	cuda_context_ewald(int ys) {
		ysize = 1;
		while (ysize < ys) {
			ysize *= 2;
		}
		CUDA_CHECK(cudaMalloc(&L, sizeof(expansion<float> )));
		CUDA_CHECK(cudaMalloc(&y, sizeof(multi_src) * ysize));
		CUDA_CHECK(cudaMallocHost(&Lp, sizeof(expansion<float> )));
		CUDA_CHECK(cudaMallocHost(&yp, sizeof(multi_src) * ysize));
		CUDA_CHECK(cudaStreamCreate(&stream));
	}
	void resize(int ys) {
		if (ys > ysize) {
			while (ysize < ys) {
				ysize *= 2;
			}
			CUDA_CHECK(cudaFree(y));
			CUDA_CHECK(cudaMalloc(&y, sizeof(multi_src) * ysize));
			CUDA_CHECK(cudaFreeHost(yp));
			CUDA_CHECK(cudaMallocHost(&yp, sizeof(multi_src) * ysize));
		}
	}
};

static std::atomic<int> lock_ewald(0);
static std::stack<cuda_context_ewald> stack_ewald;

cuda_context_ewald pop_context_ewald(int ys) {
	while (lock_ewald++ != 0) {
		lock_ewald--;
	}
	if (stack_ewald.empty()) {
		lock_ewald--;
		return cuda_context_ewald(ys);
	} else {
		auto ctx = stack_ewald.top();
		stack_ewald.pop();
		lock_ewald--;
		ctx.resize(ys);
		return ctx;
	}
}

void push_context_ewald(cuda_context_ewald ctx) {
	while (lock_ewald++ != 0) {
		lock_ewald--;
	}
	stack_ewald.push(ctx);
	lock_ewald--;
}

void gravity_CC_ewald_cuda(expansion<float> &L, const vect<position> &x, std::vector<const multi_src*> &y, bool do_phi) {

	cuda_init();

	auto ctx = pop_context_ewald(y.size());
	int k = 0;
	for (int i = 0; i < y.size(); i++) {
		ctx.yp[k++] = *y[i];
	}
	*ctx.Lp = L;
	CUDA_CHECK(cudaMemcpyAsync(ctx.y, ctx.yp, sizeof(multi_src) * y.size(), cudaMemcpyHostToDevice, ctx.stream));
	CUDA_CHECK(cudaMemcpyAsync(ctx.L, ctx.Lp, sizeof(expansion<float> ), cudaMemcpyHostToDevice, ctx.stream));

	int tb_size = (((y.size() - 1) / CCSIZE) + 1) * CCSIZE;

	/**/CC_ewald_kernel<<<dim3(tb_size/CCSIZE,1,1),dim3(CCSIZE,1,1),0,ctx.stream>>>(ctx.L, x, ctx.y, y.size(), do_phi, flop_ptr);

	CUDA_CHECK(cudaMemcpyAsync(ctx.Lp, ctx.L, sizeof(expansion<float> ), cudaMemcpyDeviceToHost, ctx.stream));
	while (cudaStreamQuery(ctx.stream) != cudaSuccess) {
		yield_to_hpx();
	}
	L = *ctx.Lp;
	push_context_ewald(std::move(ctx));
}

template<bool DO_PHI>
/**/__global__ /**/
void PP_direct_kernel(_4force *F, const vect<position> *x, const vect<position> *y, const std::pair<int, int> *yiters, int *xindex,
		int *yindex, float m, float h, double *flop_ptr) {
//	printf("sizeof(_4force) = %li\n", sizeof(_4force));

	const int iwarp = threadIdx.y;
	const int ui = blockIdx.x;
	const int l = iwarp * blockDim.x + threadIdx.x;
	const int n = threadIdx.x;
	const float Hinv = 1.0 / h;
	const float H3inv = Hinv * Hinv * Hinv;
	const float halfh = 0.5 * h;

	__shared__ vect<position>
	X[NODESIZE];
	__shared__ _4force
	G[NWARP][WARPSIZE];
	__shared__ vect<position>
	Ymem[NWARP][WARPSIZE][SYNCRATE];

	__shared__ std::uint64_t
	flop[NWARP][WARPSIZE];

	flop[iwarp][n] = 0;

	const auto yb = yindex[ui];
	const auto ye = yindex[ui + 1];
	const auto xb = xindex[ui];
	const auto xe = xindex[ui + 1];
	const auto xsize = xe - xb;
	if (l < xsize) {
		X[l] = x[xb + l];
	}
	__syncthreads();
	{
		const auto ymax = ((ye - yb - 1) / WORKSIZE + 1) * WORKSIZE + yb;
		for (int yi = yb + l; yi < ymax; yi += WORKSIZE) {
			int jb, je;
			if (yi < ye) {
				jb = yiters[yi].first;
				je = yiters[yi].second;
			}
			for (int k = 0; k < WARPSIZE; k++) {
				auto *Yptr = reinterpret_cast<float*>(Ymem[iwarp][k]);
				const int this_yi = ((yi - yb) / WARPSIZE) * WARPSIZE + k + yb;
				if (this_yi < ye) {
					const int jb = yiters[this_yi].first;
					const int je = yiters[this_yi].second;
					const int size = (je - jb) * NDIM;
					if (n < size) {
						Yptr[n] = reinterpret_cast<const float*>(y + jb)[n];
					}
				}
			}
			for (int i = xb; i < xe; i++) {
				const auto this_x = X[i - xb];
				G[iwarp][n].phi = 0.0;
				G[iwarp][n].g = vect<float>(0.0);
				if (yi < ye) {
#pragma loop unroll SYNCRATE
					for (int j0 = 0; j0 < SYNCRATE; j0++) {
						const int j = j0 + jb;
						if (j < je) {
							const vect<position> Y = Ymem[iwarp][n][j0];
							vect<float> dX;
							for (int dim = 0; dim < NDIM; dim++) {
								dX[dim] = float(int(this_x[dim]) - int(Y[dim])) * float(POS_INV);			// 3
							}
							flop[iwarp][n] += 3;
							const float r2 = dX.dot(dX);								   // 5
							const float r = sqrt(r2);// 1
							const float rinv = float(1) / max(r, halfh);// 2
							const float rinv3 = rinv * rinv * rinv;// 2
							flop[iwarp][n] += DO_PHI ? 21 : 19;
							float f, p;
							if (r > h) {
								f = rinv3;
								p = rinv;
							} else if (r > 0.5 * h) {
								const float roh = min(r * Hinv, 1.0);                         // 2
								const float roh2 = roh * roh;// 1
								const float roh3 = roh2 * roh;// 1
								f = float(-32.0 / 3.0);
								f = fma(f, roh, float(+192.0 / 5.0));// 2
								f = fma(f, roh, float(-48.0));// 2
								f = fma(f, roh, float(+64.0 / 3.0));// 2
								f = fma(f, roh3, float(-1.0 / 15.0));// 2
								f *= rinv3;// 1
								flop[iwarp][n] += 13;
								if (DO_PHI) {
									p = float(+32.0 / 15.0);
									p = fma(p, roh, float(-48.0 / 5.0));                                 // 2
									p = fma(p, roh, float(+16.0));// 2
									p = fma(p, roh, float(-32.0 / 3.0));// 2
									p = fma(p, roh2, float(+16.0 / 5.0));// 2
									p = fma(p, roh, float(-1.0 / 15.0));// 2
									p *= rinv;// 1
									flop[iwarp][n] += 11;
								}
							} else {
								const float roh = min(r * Hinv, 1.0);                           // 2
								const float roh2 = roh * roh;// 1
								f = float(+32.0);
								f = fma(f, roh, float(-192.0 / 5.0));// 2
								f = fma(f, roh2, float(+32.0 / 3.0));// 2
								f *= H3inv;// 1
								flop[iwarp][n] += 7;
								if (DO_PHI) {
									p = float(-32.0 / 5.0);
									p = fma(p, roh, float(+48.0 / 5.0));							// 2
									p = fma(p, roh2, float(-16.0 / 3.0));// 2
									p = fma(p, roh2, float(+14.0 / 5.0));// 2
									p *= Hinv;// 1
									flop[iwarp][n] += 8;
								}
							}
							const auto dXM = dX * m;								// 3
							for (int dim = 0; dim < NDIM; dim++) {
								G[iwarp][n].g[dim] -= dXM[dim] * f;    				// 6
							}
							// 13S + 2D = 15
							if( DO_PHI ) {
								G[iwarp][n].phi -= p * m;    							// 2
							}
						}
					}
				}
				for (int N = WARPSIZE / 2; N > 0; N >>= 1) {
					if (n < N) {
						G[iwarp][n].g += G[iwarp][n + N].g;
						flop[iwarp][n] += 4;
						if( DO_PHI ) {
							G[iwarp][n].phi += G[iwarp][n + N].phi;
							flop[iwarp][n] += 1;
						}
					}
				}
				if (n == 0) {
					for (int dim = 0; dim < NDIM; dim++) {
						atomicAdd(&F[i].g[dim], G[iwarp][0].g[dim]);
					}
					flop[iwarp][n] += 4;
					if( DO_PHI ) {
						atomicAdd(&F[i].phi, G[iwarp][0].phi);
						flop[iwarp][n] += 1;
					}
				}
			}
		}
	}
	for (int N = WARPSIZE / 2; N > 0; N >>= 1) {
		if (n < N) {
			flop[iwarp][n] += flop[iwarp][n + N];
		}
	}
	if (n == 0) {
		atomicAdd(flop_ptr, flop[iwarp][0]);
	}
}

__global__
void PC_direct_kernel(_4force *F, const vect<position> *x, const multi_src *z, int *xindex, int *zindex, bool do_phi, double *flop_ptr) {
//	printf("sizeof(_4force) = %li\n", sizeof(_4force));

	const int iwarp = threadIdx.y;
	const int ui = blockIdx.x;
	const int l = iwarp * blockDim.x + threadIdx.x;
	const int n = threadIdx.x;

	__shared__ vect<position>
	X[NODESIZE];
	__shared__ _4force
	G[PCNWARP][WARPSIZE];
	__shared__ std::uint64_t
	flop[NWARP][WARPSIZE];

	flop[iwarp][n] = 0;

	const auto xb = xindex[ui];
	const auto xe = xindex[ui + 1];
	const auto xsize = xe - xb;
	if (l < xsize) {
		X[l] = x[xb + l];
	}
	__syncthreads();
	int zmax = ((zindex[ui + 1] - zindex[ui] - 1) / PCWORKSIZE + 1) * PCWORKSIZE + zindex[ui];
	for (int zi = zindex[ui] + l; zi < zmax; zi += PCWORKSIZE) {
		for (int i = xb; i < xe; i++) {
			G[iwarp][n].phi = 0.0;
			G[iwarp][n].g = vect<float>(0.0);
			if (zi < zindex[ui + 1]) {
				const multipole<float> &M = z[zi].m;
				const vect<position> &Y = z[zi].x;
				vect<float> dX;
				for (int dim = 0; dim < NDIM; dim++) {
					dX[dim] = float(int(X[i - xb][dim]) - int(Y[dim])) * float(POS_INV); // 3
				}
				flop[iwarp][n] += 3;

				vect<float> g;
				float phi;
				flop[iwarp][n] += 4 + multipole_interaction(g, phi, M, dX, false, do_phi); // 516
				G[iwarp][n].g += g; // 0 / 3
				G[iwarp][n].phi += phi; // 0 / 1
			}
			for (int N = WARPSIZE / 2; N > 0; N >>= 1) {
				if (n < N) {
					G[iwarp][n].g += G[iwarp][n + N].g;
					G[iwarp][n].phi += G[iwarp][n + N].phi;
					flop[iwarp][n] += 4;
				}
			}
			if (n == 0) {
				for (int dim = 0; dim < NDIM; dim++) {
					atomicAdd(&F[i].g[dim], G[iwarp][0].g[dim]);
				}
				atomicAdd(&F[i].phi, G[iwarp][0].phi);
				flop[iwarp][0] += 4;
			}
		}
	}
	for (int N = WARPSIZE / 2; N > 0; N >>= 1) {
		if (n < N) {
			flop[iwarp][n] += flop[iwarp][n + N];
		}
	}
	if (n == 0) {
		atomicAdd(flop_ptr, flop[iwarp][0]);
	}
}

struct cuda_context {
	int xsize, ysize, zsize, isize, ypsize;
	cudaStream_t stream;
	_4force *f;
	vect<position> *x;
	std::pair<int, int> *y;
	vect<position> *ypos;
	multi_src *z;
	int *xi;
	int *yi;
	int *zi;
	_4force *fp;
	vect<position> *xp;
	multi_src *zp;
	std::pair<int, int> *yp;
	int *xip;
	int *yip;
	int *zip;
	cuda_context(int xs, int ys, int zs, int is, int yps) {
		xsize = 1;
		ysize = 1;
		zsize = 1;
		isize = 1;
		ypsize = 1;
		while (xsize < xs) {
			xsize *= 2;
		}
		while (zsize < zs) {
			zsize *= 2;
		}
		while (ysize < ys) {
			ysize *= 2;
		}
		while (ypsize < yps) {
			ypsize *= 2;
		}
		while (isize < is) {
			isize *= 2;
		}
		CUDA_CHECK(cudaMalloc(&f, sizeof(_4force) * xsize));
		CUDA_CHECK(cudaMalloc(&x, sizeof(vect<position> ) * xsize));
		CUDA_CHECK(cudaMalloc(&y, sizeof(std::pair<int, int>) * ysize));
		CUDA_CHECK(cudaMalloc(&ypos, sizeof(vect<position> ) * ypsize));
		CUDA_CHECK(cudaMalloc(&z, sizeof(multi_src) * zsize));
		CUDA_CHECK(cudaMalloc(&xi, sizeof(int) * isize));
		CUDA_CHECK(cudaMalloc(&yi, sizeof(int) * isize));
		CUDA_CHECK(cudaMalloc(&zi, sizeof(int) * isize));
		CUDA_CHECK(cudaMallocHost(&fp, sizeof(_4force) * xsize));
		CUDA_CHECK(cudaMallocHost(&xp, sizeof(vect<position> ) * xsize));
		CUDA_CHECK(cudaMallocHost(&yp, sizeof(std::pair<int, int>) * ysize));
		CUDA_CHECK(cudaMallocHost(&zp, sizeof(multi_src) * zsize));
		CUDA_CHECK(cudaMallocHost(&xip, sizeof(int) * isize));
		CUDA_CHECK(cudaMallocHost(&yip, sizeof(int) * isize));
		CUDA_CHECK(cudaMallocHost(&zip, sizeof(int) * isize));
		CUDA_CHECK(cudaStreamCreate(&stream));
	}
	void resize(int xs, int ys, int zs, int is, int yps) {
		if (yps > ypsize) {
			while (ypsize < yps) {
				ypsize *= 2;
			}
			CUDA_CHECK(cudaFree(ypos));
			CUDA_CHECK(cudaMalloc(&ypos, sizeof(vect<position> ) * ypsize));
		}
		if (xs > xsize) {
			while (xsize < xs) {
				xsize *= 2;
			}
			CUDA_CHECK(cudaFree(x));
			CUDA_CHECK(cudaFree(f));
			CUDA_CHECK(cudaMalloc(&f, sizeof(_4force) * xsize));
			CUDA_CHECK(cudaMalloc(&x, sizeof(vect<position> ) * xsize));
			CUDA_CHECK(cudaFreeHost(xp));
			CUDA_CHECK(cudaFreeHost(fp));
			CUDA_CHECK(cudaMallocHost(&fp, sizeof(_4force) * xsize));
			CUDA_CHECK(cudaMallocHost(&xp, sizeof(vect<position> ) * xsize));
		}
		if (ys > ysize) {
			while (ysize < ys) {
				ysize *= 2;
			}
			CUDA_CHECK(cudaFree(y));
			CUDA_CHECK(cudaMalloc(&y, sizeof(std::pair<int, int>) * ysize));
			CUDA_CHECK(cudaFreeHost(yp));
			CUDA_CHECK(cudaMallocHost(&yp, sizeof(std::pair<int, int>) * ysize));
		}
		if (zs > zsize) {
			while (zsize < zs) {
				zsize *= 2;
			}
			CUDA_CHECK(cudaFree(z));
			CUDA_CHECK(cudaMalloc(&z, sizeof(multi_src) * zsize));
			CUDA_CHECK(cudaFreeHost(zp));
			CUDA_CHECK(cudaMallocHost(&zp, sizeof(multi_src) * zsize));
		}
		if (is > isize) {
			while (isize < is) {
				isize *= 2;
			}
			CUDA_CHECK(cudaFree(xi));
			CUDA_CHECK(cudaFree(yi));
			CUDA_CHECK(cudaFree(zi));
			CUDA_CHECK(cudaMalloc(&xi, sizeof(int) * isize));
			CUDA_CHECK(cudaMalloc(&yi, sizeof(int) * isize));
			CUDA_CHECK(cudaMalloc(&zi, sizeof(int) * isize));
			CUDA_CHECK(cudaFreeHost(xip));
			CUDA_CHECK(cudaFreeHost(yip));
			CUDA_CHECK(cudaFreeHost(zip));
			CUDA_CHECK(cudaMallocHost(&xip, sizeof(int) * isize));
			CUDA_CHECK(cudaMallocHost(&yip, sizeof(int) * isize));
			CUDA_CHECK(cudaMallocHost(&zip, sizeof(int) * isize));
		}
	}
};

static std::atomic<int> lock(0);
static std::stack<cuda_context> stack;

cuda_context pop_context(int xs, int ys, int zs, int is, int yps) {
	while (lock++ != 0) {
		lock--;
	}
	if (stack.empty()) {
		lock--;
		return cuda_context(xs, ys, zs, is, yps);
	} else {
		auto ctx = stack.top();
		stack.pop();
		lock--;
		ctx.resize(xs, ys, zs, is, yps);
		return ctx;
	}
}

void push_context(cuda_context ctx) {
	while (lock++ != 0) {
		lock--;
	}
	stack.push(ctx);
	lock--;
}

struct pinned_context {
	pinned_vector<int> xindex;
	pinned_vector<int> yindex;
	pinned_vector<_4force> f;
	pinned_vector<vect<position>> x;
	pinned_vector<std::pair<int, int>> y;
	pinned_vector<multi_src> z;
	pinned_vector<int> zindex;
	pinned_context() = default;
	pinned_context(const pinned_context&) = delete;
	pinned_context(pinned_context&&) = default;
	pinned_context& operator=(const pinned_context&) = delete;
	pinned_context& operator=(pinned_context&&) = default;
};

std::stack<pinned_context> pinned_stack;
std::atomic<int> pinned_mtx(0);

pinned_context pop_pinned_context() {
	while (pinned_mtx++ != 0) {
		pinned_mtx--;
	}
	pinned_context ctx;
	if (!pinned_stack.empty()) {
		ctx = std::move(pinned_stack.top());
		pinned_stack.pop();
	}
	pinned_mtx--;
	return std::move(ctx);
}

void push_pinned_context(pinned_context &&ctx) {
	while (pinned_mtx++ != 0) {
		pinned_mtx--;
	}
	pinned_stack.push(std::move(ctx));
	pinned_mtx--;
}

void gravity_PP_direct_cuda(std::vector<cuda_work_unit> &&units, const pinned_vector<vect<position>> &ydata, bool do_phi) {
	static const float m = opts.particle_mass;
	cuda_init();
	std::uint64_t interactions = 0;
	{
		auto pctx = pop_pinned_context();
		auto &xindex = pctx.xindex;
		auto &yindex = pctx.yindex;
		auto &f = pctx.f;
		auto &x = pctx.x;
		auto &y = pctx.y;
		xindex.resize(0);
		yindex.resize(0);
		f.resize(0);
		x.resize(0);
		y.resize(0);
		int xi = 0;
		int yi = 0;
		for (const auto &unit : units) {
			xindex.push_back(xi);
			yindex.push_back(yi);
			xi += unit.xptr->size();
			yi += unit.yiters.size();
			for (const auto &this_f : *unit.fptr) {
				f.push_back(this_f);
			}
			for (const auto &this_x : *unit.xptr) {
				x.push_back(this_x);
			}
			for (int j = 0; j < unit.yiters.size(); j++) {
				std::pair<int, int> iter = unit.yiters[j];
				interactions += unit.xptr->size() * (iter.second - iter.first);
				y.push_back(iter);
			}
		}
		xindex.push_back(xi);
		yindex.push_back(yi);
		const auto fbytes = sizeof(_4force) * f.size();
		const auto xbytes = sizeof(vect<position> ) * x.size();
		const auto ybytes = sizeof(std::pair<int, int>) * y.size();
		const auto ypbytes = sizeof(vect<position> ) * ydata.size();
		const auto xibytes = sizeof(int) * xindex.size();
		const auto yibytes = sizeof(int) * yindex.size();

		auto ctx = pop_context(x.size(), y.size(), 0, xindex.size(), ydata.size());
		memcpy(ctx.fp, f.data(), fbytes);
		memcpy(ctx.xp, x.data(), xbytes);
		memcpy(ctx.yp, y.data(), ybytes);
		memcpy(ctx.xip, xindex.data(), xibytes);
		memcpy(ctx.yip, yindex.data(), yibytes);
		CUDA_CHECK(cudaMemcpyAsync(ctx.f, ctx.fp, fbytes, cudaMemcpyHostToDevice, ctx.stream));
		CUDA_CHECK(cudaMemcpyAsync(ctx.y, y.data(), ybytes, cudaMemcpyHostToDevice, ctx.stream));
		CUDA_CHECK(cudaMemcpyAsync(ctx.ypos, ydata.data(), ypbytes, cudaMemcpyHostToDevice, ctx.stream));
		CUDA_CHECK(cudaMemcpyAsync(ctx.x, ctx.xp, xbytes, cudaMemcpyHostToDevice, ctx.stream));
		CUDA_CHECK(cudaMemcpyAsync(ctx.yi, yindex.data(), yibytes, cudaMemcpyHostToDevice, ctx.stream));
		CUDA_CHECK(cudaMemcpyAsync(ctx.xi, xindex.data(), xibytes, cudaMemcpyHostToDevice, ctx.stream));
		if (do_phi) {
		PP_direct_kernel<true><<<dim3(units.size(),1,1),dim3(WARPSIZE,NWARP,1),0,ctx.stream>>>(ctx.f,ctx.x,ctx.ypos, ctx.y,ctx.xi,ctx.yi, m, opts.soft_len, flop_ptr);
	} else {
	PP_direct_kernel<false><<<dim3(units.size(),1,1),dim3(WARPSIZE,NWARP,1),0,ctx.stream>>>(ctx.f,ctx.x,ctx.ypos, ctx.y,ctx.xi,ctx.yi, m, opts.soft_len, flop_ptr);
}

CUDA_CHECK(cudaMemcpyAsync(f.data(), ctx.f, fbytes, cudaMemcpyDeviceToHost, ctx.stream));
while (cudaStreamQuery(ctx.stream) != cudaSuccess) {
	yield_to_hpx();
}
int k = 0;
for (const auto &unit : units) {
	for (auto &this_f : *unit.fptr) {
		this_f = f[k];
		k++;
	}
}
push_context(std::move(ctx));
push_pinned_context(std::move(pctx));

}
{
auto pctx = pop_pinned_context();
auto &xindex = pctx.xindex;
auto &z = pctx.z;
auto &zindex = pctx.zindex;
auto &f = pctx.f;
auto &x = pctx.x;
xindex.resize(0);
zindex.resize(0);
f.resize(0);
x.resize(0);
z.resize(0);

int xi = 0;
int zi = 0;
int size = 0;
std::uint64_t interactions = 0;
for (const auto &unit : units) {
	if (unit.z.size()) {
		xindex.push_back(xi);
		zindex.push_back(zi);
		xi += unit.xptr->size();
		zi += unit.z.size();
		for (const auto &this_f : *unit.fptr) {
			f.push_back(this_f);
		}
		for (const auto &this_x : *unit.xptr) {
			x.push_back(this_x);
		}
		for (int j = 0; j < unit.z.size(); j++) {
			z.push_back(*unit.z[j]);
		}
		size++;
	}
}
xindex.push_back(xi);
zindex.push_back(zi);
if (z.size()) {
	const auto fbytes = sizeof(_4force) * f.size();
	const auto xbytes = sizeof(vect<position> ) * x.size();
	const auto zbytes = sizeof(multi_src) * z.size();
	const auto xibytes = sizeof(int) * xindex.size();
	const auto zibytes = sizeof(int) * zindex.size();

	auto ctx = pop_context(x.size(), 0, z.size(), zindex.size(), 0);
	CUDA_CHECK(cudaMemcpyAsync(ctx.f, f.data(), fbytes, cudaMemcpyHostToDevice, ctx.stream));
//		printf( "%li %lli %lli\n", zbytes, ctx.z, ctx.zp);
	CUDA_CHECK(cudaMemcpyAsync(ctx.z, z.data(), zbytes, cudaMemcpyHostToDevice, ctx.stream));
	CUDA_CHECK(cudaMemcpyAsync(ctx.x, x.data(), xbytes, cudaMemcpyHostToDevice, ctx.stream));
	CUDA_CHECK(cudaMemcpyAsync(ctx.xi, xindex.data(), xibytes, cudaMemcpyHostToDevice, ctx.stream));
	CUDA_CHECK(cudaMemcpyAsync(ctx.zi, zindex.data(), zibytes, cudaMemcpyHostToDevice, ctx.stream));

	/**/PC_direct_kernel<<<dim3(size,1,1),dim3(WARPSIZE,PCNWARP,1),0,ctx.stream>>>(ctx.f,ctx.x,ctx.z,ctx.xi,ctx.zi, do_phi, flop_ptr);

			CUDA_CHECK(cudaMemcpyAsync(f.data(), ctx.f, fbytes, cudaMemcpyDeviceToHost, ctx.stream));
	while (cudaStreamQuery(ctx.stream) != cudaSuccess) {
		yield_to_hpx();
	}
	int k = 0;
	for (const auto &unit : units) {
		if (unit.z.size()) {
			for (auto &this_f : *unit.fptr) {
				this_f = f[k];
				k++;
			}
		}
	}
	push_context(ctx);
}
push_pinned_context(std::move(pctx));
}
}

