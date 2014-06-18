/* This is a CUDA implementation of Jermann Quadrini 2013 AER  
 * We simulate a equilibrium that collateral constraint is always binding to check the accuracy of 
 * their linearization approach. Hopeully we can ind something that they missed. A main suspect 
 * is the asymmetry of policy functions.
 */

#define nk 2560
#define nz 9
#define nxxi 9
#define nm1 5120 
#define tol 1e-6
#define maxiter 2500
#define kwidth 1.2
#define mkwidth 20.0 

/* Includes, system */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// Includes, Thrust
#include <thrust/functional.h>
#include <thrust/for_each.h>
#include <thrust/sort.h>
#include <thrust/extrema.h>
#include <thrust/tuple.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sequence.h>
#include <thrust/device_ptr.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/zip_iterator.h>

/* Includes, cuda */
// #include <cuda_runtime.h>
#include <cublas_v2.h>
// #include <helper_functions.h>
#include "cuda_helpers.h"

// Includes, C++ codes
#include "cppcode.h"

// Includes, model specific things
#include "model.h"

using namespace std;
using namespace thrust;

void guess_vfi(const host_vector<double> K, const host_vector<double> Z, const host_vector<double> XXI, host_vector<double> & V1_low, host_vector<double> & V1_high, para p, double factor) {
	// Try results from vfi iteration
	host_vector<double> K_vfi;
	host_vector<double> Z_vfi;
	host_vector<double> XXI_vfi;
	host_vector<double> M_vfi;
	load_vec(K_vfi,"./vfi_results/Kgrid.csv");
	load_vec(Z_vfi,"./vfi_results/Zgrid.csv");
	load_vec(XXI_vfi,"./vfi_results/XXIgrid.csv");
	load_vec(M_vfi,"./vfi_results/mopt.csv");

	// Create guesses.
	for (int i_k=0; i_k<K.size(); i_k++) {
		int i_k_vfi = fit2grid(K[i_k],K_vfi);
		for (int i_z = 0; i_z < Z.size(); i_z++) {
			int i_z_vfi = fit2grid(Z[i_z],Z_vfi);
			for (int i_xxi = 0; i_xxi < XXI.size(); i_xxi++) {
				int i_xxi_vfi = fit2grid(XXI[i_xxi],XXI_vfi);
				double temp = M_vfi[i_k_vfi+i_z_vfi*K_vfi.size()+i_xxi_vfi*K_vfi.size()*Z_vfi.size()];
				V1_low[i_k+nk*i_z+nk*nz*i_xxi] = 1/factor*temp;
				V1_high[i_k+nk*i_z+nk*nz*i_xxi] = factor*temp;
			};
		};
	};
};

void guess_linear(const host_vector<double> K, const host_vector<double> Z, const host_vector<double> XXI, host_vector<double> & V1_low, host_vector<double> & V1_high, para p, double factor_low, double factor_high) {
	// Initialize matrices
	int n = 9; int n_jump = 8; int n_shock = 2;
	host_vector<double> A(n*n,0); 
	host_vector<double> B(n*n,0); 
	host_vector<double> C(n*n_shock,0); 
	host_vector<double> rrho(n_shock*n_shock,0);
   	host_vector<double> Pphi(n*(n-n_jump+n_shock),0);

	// Fill in matrices.
	linearizedmodel(A.data(),B.data(),C.data(),rrho.data(),n,n_shock,p);

	// Call linear solver
	linearQZ(A.data(),B.data(),C.data(),rrho.data(),n,n_jump,n_shock,Pphi.data());

	// Create guesses.
	for (int i_k=0; i_k<nk; i_k++) {
		for (int i_z = 0; i_z < nz; i_z++) {
			for (int i_xxi = 0; i_xxi < nxxi; i_xxi++) {
				double temp = p.mkss+Pphi[8+0*9]*(K[i_k]-p.kss) + Pphi[8+1*9]*(log(Z[i_z])-log(p.zbar))+ Pphi[8+2*9]*(log(XXI[i_xxi])-log(p.xxibar));
				V1_low[i_k+nk*i_z+nk*nz*i_xxi] = factor_low*temp;
				V1_high[i_k+nk*i_z+nk*nz*i_xxi] = factor_high*temp;
			};
		};
	};
};

__host__ __device__
bool eureka(state s, shadow m, control & u1, control & u2, para p, int i_z, int i_xxi, double* EM1_low, double* EM1_high, double* K) {
	double interp_low, interp_high;
	int i_kplus;

	// Case 1: Binding 
	u1.compute(s,m,p,1);

	// A series of tests whether it make senses
	if (u1.c <= 0) {
		goto case2;
	};
	if (u1.kplus <= 0) {
		goto case2;
	};
	// if (u.kplus < k[0]) {
	// 	goto case2;
	// };
	// if (u.kplus > k[nk-1]) {
	// 	goto case2;
	// };
	if (u1.mmu < 0) {
		// LM shouldn't be negative
		goto case2;
	};
	if ( (u1.n <= 0) || (u1.n >= 1) ) {
		// Hours out of bound.
		goto case2;
	};

	i_kplus = fit2grid(u1.kplus,nk,K);
	interp_low = EM1_low[i_kplus+i_z*nk+i_xxi*nz*nk]+(u1.kplus-K[i_kplus])*(EM1_low[i_kplus+1+i_z*nk+i_xxi*nz*nk]-EM1_low[i_kplus+i_z*nk+i_xxi*nz*nk])/(K[i_kplus+1]-K[i_kplus]);
	interp_high = EM1_high[i_kplus+i_z*nk+i_xxi*nz*nk]+(u1.kplus-K[i_kplus])*(EM1_high[i_kplus+1+i_z*nk+i_xxi*nz*nk]-EM1_high[i_kplus+i_z*nk+i_xxi*nz*nk])/(K[i_kplus+1]-K[i_kplus]);
	if ( (u1.lhs1 > p.bbeta*interp_high) || (p.bbeta*interp_low > u1.lhs1) ) {
		// Euler equation 1 fails
		goto case2;
	};

	// Found a sensible m, break the loop
	return true;

case2: // Not Binding
	u2.compute(s,m,p,0);
	// A series of tests whether it make senses
	if (u2.c <= 0) {
		return false;
	};
	if (u2.kplus <= 0) {
		return false;
	};
	// if (u2.kplus < K[0]) {
	// 	return false;
	// };
	// if (u2.kplus > K[nk-1]) {
	// 	return false;
	// };
	if (s.xxi*u2.kplus <= u2.Y) {
		// Financial constraint shouldn't be binding 
		return false;
	};
	if ( (u2.n <= 0) || (u2.n >= 1) ) {
		// Hours out of bound.
		return false;
	};

	i_kplus = fit2grid(u2.kplus,nk,K);
	interp_low = EM1_low[i_kplus+i_z*nk+i_xxi*nz*nk]+(u2.kplus-K[i_kplus])*(EM1_low[i_kplus+1+i_z*nk+i_xxi*nz*nk]-EM1_low[i_kplus+i_z*nk+i_xxi*nz*nk])/(K[i_kplus+1]-K[i_kplus]);
	interp_high = EM1_high[i_kplus+i_z*nk+i_xxi*nz*nk]+(u2.kplus-K[i_kplus])*(EM1_high[i_kplus+1+i_z*nk+i_xxi*nz*nk]-EM1_high[i_kplus+i_z*nk+i_xxi*nz*nk])/(K[i_kplus+1]-K[i_kplus]);
	if ( (u2.lhs1 > p.bbeta*interp_high) || (p.bbeta*interp_low > u2.lhs1) ) {
		// Euler equation 1 fails
		return false;
	};

	// A sensible m found, we escape the loop
	return true;
};

// This functor yields RHS for each (k', k, z). Follwing examples in Thrust doc
struct shrink 
{
	// Data member
	double *K, *Z, *XXI;
	double *V1_low;
	double *V1_high;
	double *Vplus1_low;
	double *Vplus1_high;
	double *EM1_low;
	double *EM1_high;
	double *flag;
	para p;

	// Construct this object, create util from _util, etc.
	__host__ __device__
	shrink(double* K_ptr, double* Z_ptr, double* XXI_ptr,
	double* V1_low_ptr,
	double* V1_high_ptr,
	double* Vplus1_low_ptr,
	double* Vplus1_high_ptr,
	double* EM1_low_ptr,
	double* EM1_high_ptr,
	double* flag_ptr,
	para _p)
	{
		K = K_ptr; Z = Z_ptr; XXI = XXI_ptr;
		V1_low = V1_low_ptr;
		V1_high = V1_high_ptr;
		Vplus1_low = Vplus1_low_ptr;
		Vplus1_high = Vplus1_high_ptr;
		EM1_low = EM1_low_ptr;
		EM1_high = EM1_high_ptr;
		flag = flag_ptr;
		p = _p;
	};

	__host__ __device__
	void operator()(int index) {
		// Perform ind2sub
		int subs[3];
		int size_vec[3];
		size_vec[0] = nk;
		size_vec[1] = nz;
		size_vec[2] = nxxi;
		ind2sub(3,size_vec,index,subs);
		int i_k = subs[0];
		int i_z = subs[1];
		int i_xxi = subs[2];

		// Find the "box" or "hypercube" that described m's range. Fancy word.
		double m1min = V1_low[index]; double m1max = V1_high[index];
		double m1min_old = m1min;
		double m1max_old = m1max;
		double step1 = (m1max_old-m1min_old)/double(nm1-1);
		double tempflag = 0.0;

		// Find and construct state and control, otherwise they won't update in the for loop
		double k =K[i_k]; double z=Z[i_z]; double xxi=XXI[i_xxi];
		state s(k,z,xxi,p);
		control u1, u2;

		// Initial search to find the min m
		for (int i_m1min = 0; i_m1min < nm1; i_m1min++) {
			// Construct state and find control variables
			double m1 = m1min+double(i_m1min)*step1;
			if (eureka(s,shadow(m1),u1,u2,p,i_z,i_xxi,EM1_low,EM1_high,K)) {
				m1min = m1;
				tempflag++;
				break;
			};
		};

		// Initial search to find the max m
		for (int i_m1max = 0; i_m1max < nm1; i_m1max++) {
			// Construct state and find control variables
			double m1 = m1max - double(i_m1max)*step1;
			if (eureka(s,shadow(m1),u1,u2,p,i_z,i_xxi,EM1_low,EM1_high,K)) {
				m1max = m1;
				tempflag++;
				break;
			};
		};

		// Update Vs
		flag[index] = double(tempflag)/double(nm1);
		Vplus1_high[index] = m1max + step1*int((m1max<m1max_old));
		Vplus1_low[index] = m1min - step1*int((m1min>m1min_old));
	};
};	

// This functor calculates the error
struct myMinus {
	// Tuple is (V1low,Vplus1low,V1high,Vplus1high,...)
	template <typename Tuple>
	__host__ __device__
	double operator()(Tuple t)
	{
		return max( abs(get<0>(t)-get<1>(t)),abs(get<2>(t)-get<3>(t)) );
	}
};

// This functor calculates the distance 
struct myDist {
	// Tuple is (V1low,Vplus1low,V1high,Vplus1high,...)
	template <typename Tuple>
	__host__ __device__
	double operator()(Tuple t)
	{
		return abs(get<0>(t)-get<1>(t));
	}
};

int main(int argc, char ** argv)
{
	// Initialize Parameters
	para p;

	// Set Model Parameters
	p.bbeta = 0.9825;
	p.ddelta = 0.025;
	p.ttheta = 0.36;
	p.kkappa = 0.1460;
	p.ttau = 0.3500;
	p.xxibar = 0.11;
	p.zbar = 1.0;
	p.rrhozz = 0.9457;
	p.rrhoxxiz = 0.0321;
	p.rrhozxxi =-0.0091;
	p.rrhoxxixxi = 0.9703;
	p.var_epsz = 0.0045*0.0045;
	p.var_epsxxi = 0.0098*0.0098;
	p.complete(); // complete all implied p. find S-S

	cout << setprecision(16) << "kss: " << p.kss << endl;
	cout << setprecision(16) << "zss: " << p.zbar << endl;
	cout << setprecision(16) << "xxiss: " <<p.xxibar << endl;
	cout << setprecision(16) << "mkss: " << p.mkss << endl;
	cout << setprecision(16) << "dss: " << p.dss << endl;
	cout << setprecision(16) << "css: " << p.css << endl;
	cout << setprecision(16) << "nss: " << p.nss << endl;
	cout << setprecision(16) << "wss: " << p.wss << endl;
	cout << setprecision(16) << "mmuss: " << p.mmuss << endl;
	cout << setprecision(16) << "aalpha: " << p.aalpha << endl;
	cout << setprecision(16) << "tol: " << tol << endl;

	// Select Device
	// int num_devices;
	// cudaGetDeviceCount(&num_devices);
	if (argc > 1) {
		int gpu = atoi(argv[1]);
		cudaSetDevice(gpu);
	};
	// Only for cuBLAS
	const double alpha = 1.0;
	const double beta = 0.0;

	// Create all STATE, SHOCK grids here
	host_vector<double> h_K(nk); 
	host_vector<double> h_Z(nz);
	host_vector<double> h_XXI(nxxi);

	host_vector<double> h_V1_low(nk*nz*nxxi, 1/mkwidth*p.mkss);
	host_vector<double> h_V1_high(nk*nz*nxxi,mkwidth*p.mkss);
	host_vector<double> h_Vplus1_low(nk*nz*nxxi,1/mkwidth*p.mkss);
	host_vector<double> h_Vplus1_high(nk*nz*nxxi,mkwidth*p.mkss);

	host_vector<double> h_EM1_low(nk*nz*nxxi,0.0);
	host_vector<double> h_EM1_high(nk*nz*nxxi,0.0);

	host_vector<double> h_P(nz*nxxi*nz*nxxi, 0);
	host_vector<double> h_flag(nk*nz*nxxi, 0); 

	// host_vector<double> h_c(nk*nz*nxxi*nm1);
	// host_vector<double> h_n(nk*nz*nxxi*nm1);
	// host_vector<double> h_kplus(nk*nz*nxxi*nm1);
	// host_vector<double> h_mmu(nk*nz*nxxi*nm1);
	// host_vector<double> h_lhs(nk*nz*nxxi*nm1);
	// host_vector<double> h_rhs_low(nk*nz*nxxi*nm1);
	// host_vector<double> h_rhs_high(nk*nz*nxxi*nm1);
	
	// Create capital grid
	double minK = 1/kwidth*p.kss;
	double maxK = kwidth*p.kss;
	linspace(minK,maxK,nk,raw_pointer_cast(h_K.data()));

	// Create shocks grids
	host_vector<double> h_shockgrids(2*nz);
	double* h_shockgrids_ptr = raw_pointer_cast(h_shockgrids.data());
	double* h_P_ptr = raw_pointer_cast(h_P.data());
	gridgen_fptr linspace_fptr = &linspace; // select linspace as grid gen
	tauchen_vec(2,nz,5,p.A,p.Ssigma_e,h_shockgrids_ptr,h_P_ptr,linspace_fptr);
	for (int i_shock = 0; i_shock < nz; i_shock++) {
		h_Z[i_shock] = p.zbar*exp(h_shockgrids[i_shock+0*nz]);
		h_XXI[i_shock] = p.xxibar*exp(h_shockgrids[i_shock+1*nz]);
	};

	// Obtain initial guess from linear solution
	guess_vfi(h_K, h_Z, h_XXI, h_V1_low, h_V1_high, p, 1.05) ;

	// Copy to the device
	device_vector<double> d_K = h_K;
	device_vector<double> d_Z = h_Z;
	device_vector<double> d_XXI = h_XXI;

	device_vector<double> d_V1_low = h_V1_low;
	device_vector<double> d_V1_high = h_V1_high;

	device_vector<double> d_Vplus1_low = h_Vplus1_low;
	device_vector<double> d_Vplus1_high = h_Vplus1_high;

	device_vector<double> d_EM1_low = h_EM1_low;
	device_vector<double> d_EM1_high = h_EM1_high;

	device_vector<double> d_P = h_P;
	device_vector<double> d_flag = h_flag;

	// Obtain device pointers to be used by cuBLAS
	double* d_K_ptr = raw_pointer_cast(d_K.data());
	double* d_Z_ptr = raw_pointer_cast(d_Z.data());
	double* d_XXI_ptr = raw_pointer_cast(d_XXI.data());

	double* d_V1_low_ptr = raw_pointer_cast(d_V1_low.data());
	double* d_V1_high_ptr = raw_pointer_cast(d_V1_high.data());

	double* d_Vplus1_low_ptr = raw_pointer_cast(d_Vplus1_low.data());
	double* d_Vplus1_high_ptr = raw_pointer_cast(d_Vplus1_high.data());

	double* d_EM1_low_ptr = raw_pointer_cast(d_EM1_low.data());
	double* d_EM1_high_ptr = raw_pointer_cast(d_EM1_high.data());

	double* d_P_ptr = raw_pointer_cast(d_P.data());
	double* d_flag_ptr = raw_pointer_cast(d_flag.data());

	// Firstly a virtual index array from 0 to nk*nk*nz
	counting_iterator<int> begin(0);
	counting_iterator<int> end(nk*nz*nxxi);

	// Start Timer
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,NULL);
	
	// Step.1 Has to start with this command to create a handle
	cublasHandle_t handle;

	// Step.2 Initialize a cuBLAS context using Create function,
	// and has to be destroyed later
	cublasCreate(&handle);
	
	double diff = 10; double dist = 100; int iter = 0;
	while ((diff>tol)&&(iter<maxiter)){
		// Find EMs for low and high 
		cublasDgemm(handle,
			CUBLAS_OP_N,  
			CUBLAS_OP_T,
			nk, nz*nxxi, nz*nxxi,
			&alpha,
			d_V1_low_ptr, 
			nk, 
			d_P_ptr,
			nz*nxxi,
			&beta,
			d_EM1_low_ptr,
			nk);
		cublasDgemm(handle,
			CUBLAS_OP_N,  
			CUBLAS_OP_T,
			nk, nz*nxxi, nz*nxxi,
			&alpha,
			d_V1_high_ptr, 
			nk, 
			d_P_ptr,
			nz*nxxi,
			&beta,
			d_EM1_high_ptr,
			nk);

		// Directly find the new Value function
		thrust::for_each(
			begin,
			end,
			shrink(d_K_ptr, d_Z_ptr, d_XXI_ptr,
				d_V1_low_ptr,
				d_V1_high_ptr,
				d_Vplus1_low_ptr,
				d_Vplus1_high_ptr,
				d_EM1_low_ptr,
				d_EM1_high_ptr,
				d_flag_ptr,
				p)
		);

		// Find error
		double diff1 = transform_reduce(
			make_zip_iterator(make_tuple(d_V1_low.begin(), d_Vplus1_low.begin(), d_V1_high.begin(),d_Vplus1_high.begin())),
			make_zip_iterator(make_tuple(d_V1_low.end()  , d_Vplus1_low.end()  , d_V1_high.end()  ,d_Vplus1_high.end())),
			myMinus(),
			0.0,
			maximum<double>()
			);

		// Find distance 
		double dist1 = transform_reduce(
			make_zip_iterator(make_tuple(d_Vplus1_low.begin(),d_Vplus1_high.begin())),
			make_zip_iterator(make_tuple(d_Vplus1_low.end()  ,d_Vplus1_high.end())),
			myDist(),
			0.0,
			maximum<double>()
			);
		diff = max(diff1,-99.0);
		dist = max(dist1,-99.0);

		cout << "diff is: "<< diff << endl;
		cout << "dist is: "<< dist << endl;
		cout << "Vplus1[117,2,4] (the spike) range is " << d_Vplus1_low[117+2*nk+4*nk*nz] << ", " << d_Vplus1_high[117+2*nk+4*nk*nz] << endl;

		// update correspondence
		d_V1_low = d_Vplus1_low; d_V1_high = d_Vplus1_high;

		cout << ++iter << endl;
		cout << "=====================" << endl;

	};

	//==========cuBLAS stuff ends=======================
	// Step.3 Destroy the handle.
	cublasDestroy(handle);

	// Stop Timer
	cudaEventRecord(stop,NULL);
	cudaEventSynchronize(stop);
	float msecTotal = 0.0;
	cudaEventElapsedTime(&msecTotal, start, stop);

	// Compute and print the performance
	float msecPerMatrixMul = msecTotal;
	cout << "Time= " << msecPerMatrixMul << " msec, iter= " << iter << endl;

	// Copy back to host and print to file
	h_V1_low = d_V1_low; h_V1_high = d_V1_high;
	h_EM1_low = d_EM1_low; h_EM1_high = d_EM1_high;
	h_flag = d_flag;
	
	// Compute and save the decision variables
	host_vector<double> h_copt(nk*nz*nxxi);
	host_vector<double> h_kopt(nk*nz*nxxi);
	host_vector<double> h_nopt(nk*nz*nxxi);
	host_vector<double> h_mmuopt(nk*nz*nxxi);
	host_vector<double> h_dopt(nk*nz*nxxi);
	host_vector<double> h_wopt(nk*nz*nxxi);
	host_vector<double> h_kk_1(nm1);
	host_vector<double> h_kk_2(nm1);
	host_vector<double> h_lhs1_1(nm1);
	host_vector<double> h_lhs1_2(nm1);
	host_vector<double> h_rhslow_1(nm1);
	host_vector<double> h_rhshigh_1(nm1);
	host_vector<double> h_rhslow_2(nm1);
	host_vector<double> h_rhshigh_2(nm1);
	host_vector<double> h_nn_1(nm1);
	host_vector<double> h_nn_2(nm1);

	for (int i_k=0; i_k<nk; i_k++) {
		for (int i_z = 0; i_z < nz; i_z++) {
			for (int i_xxi=0; i_xxi < nxxi; i_xxi++) {
				int index = i_k+i_z*nk+i_xxi*nk*nz;
				double m1 = (h_V1_high[index]+h_V1_high[index])/2;
				double k = h_K[i_k];
				double z=h_Z[i_z]; double xxi=h_XXI[i_xxi];
				control u;
				state s(k,z,xxi,p);

				// Try not binding first
				u.compute(s,shadow(m1),p,0);
				if (
						(s.xxi*u.kplus > u.Y) &&
						(u.c > 0) && 
						(u.kplus > 0) &&
						(u.n > 0) &&
						(u.n < 1) 
				   )
				{
					h_copt[index] = u.c;
					h_kopt[index] = u.kplus;
					h_nopt[index] = u.n;
					h_mmuopt[index] = u.mmu;
					h_dopt[index] = u.d;
					h_wopt[index] = u.w;
				} else {
					u.compute(s,shadow(m1),p,1);
					h_copt[index] = u.c;
					h_kopt[index] = u.kplus;
					h_nopt[index] = u.n;
					h_mmuopt[index] = u.mmu;
					h_dopt[index] = u.d;
					h_wopt[index] = u.w;
				};

			};
		};
	};
	
	save_vec(h_K,nk,"./adrian_results/Kgrid.csv");
	save_vec(h_Z,"./adrian_results/Zgrid.csv");
	save_vec(h_XXI,"./adrian_results/XXIgrid.csv");
	save_vec(h_P,"./adrian_results/P.csv");
	save_vec(h_V1_low,"./adrian_results/V1_low_guess.csv");
	save_vec(h_V1_high,"./adrian_results/V1_high_guess.csv");
	save_vec(h_V1_low,"./adrian_results/V1_low.csv");
	save_vec(h_V1_high,"./adrian_results/V1_high.csv");
	save_vec(h_flag,"./adrian_results/flag.csv");
	save_vec(h_copt,"./adrian_results/copt.csv");
	save_vec(h_kopt,"./adrian_results/kopt.csv");
	save_vec(h_nopt,"./adrian_results/nopt.csv");
	save_vec(h_mmuopt,"./adrian_results/mmuopt.csv");
	save_vec(h_dopt,"./adrian_results/dopt.csv");
	save_vec(h_wopt,"./adrian_results/wopt.csv");
	save_vec(h_kk_1,"./adrian_results/kk_1.csv");
	save_vec(h_kk_2,"./adrian_results/kk_2.csv");
	save_vec(h_nn_1,"./adrian_results/nn_1.csv");
	save_vec(h_nn_2,"./adrian_results/nn_2.csv");
	save_vec(h_lhs1_1,"./adrian_results/lhs1_1.csv");
	save_vec(h_lhs1_2,"./adrian_results/lhs1_2.csv");
	save_vec(h_rhslow_1,"./adrian_results/rhslow_1.csv");
	save_vec(h_rhshigh_1,"./adrian_results/rhshigh_1.csv");
	save_vec(h_rhslow_2,"./adrian_results/rhslow_2.csv");
	save_vec(h_rhshigh_2,"./adrian_results/rhshigh_2.csv");

	// Export parameters to MATLAB
	p.exportmatlab("./MATLAB/mypara.m");

	return 0;
}
