#define nk 25000
#define nm1 2560
#define nz 23
#define nxxi 23
#define tol 1e-7
#define maxiter 2500
#define kwidth 1.5

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

// Includes, cuda 
#include <cublas_v2.h>
#include "cuda_helpers.h"

// Includes, C++ codes
#include "cppcode.h"

// Includes, model
#include "model.h"

using namespace std;
using namespace thrust;

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

/// This function finds the value of RHS given k', k, z, xxi
__host__ __device__
double rhsvalue (state s, double kplus, int i_z, int i_xxi, int i_kplus, double* EV, para p, int & binding) {
	double temp1 = -9999999;
	double temp2 = -9999999;
	control u1, u2;
	state splus(kplus,1,1,p);
	u1.compute(s,splus,p,1);
	u2.compute(s,splus,p,0);
	if (
			// (u1.mmu >= 0) &&
			(u1.c > 0 ) &&
			(u1.n > 0 ) &&
			(u1.n < 1 )
	   )
	{
		temp1 = log(u1.c) + p.aalpha*log(1-u1.n) + p.bbeta*EV[i_kplus+i_z*nk+i_xxi*nk*nz];
	};
	if (
			(s.xxi*splus.k > u2.Y) &&
			(u2.c > 0 ) &&
			(u2.n > 0 ) &&
			(u2.n < 1 )
	   )
	{
		temp2 = log(u2.c) + p.aalpha*log(1-u2.n) + p.bbeta*EV[i_kplus+i_z*nk+i_xxi*nk*nz];
	};
	if (temp1 > temp2) {
		binding = 1;
		return temp1;
	} else {
		binding = 0;
		return temp2;
	};
};
// This find the max using binary search and assumes concavity
__host__ __device__
void concavemax(double k, double z, double xxi, const int left_ind, const int right_ind, const int i_k,const int i_z, const int i_xxi, 
		double* K, double* EV, double* Vplus, int* koptind, para p) {
	int index = i_k + i_z*nk + i_xxi*nk*nz;
	int trash;

	if (right_ind-left_ind==1) {
		double left_value, right_value;
		left_value = rhsvalue(state(k,z,xxi,p),K[left_ind],i_z,i_xxi,left_ind, EV, p, trash);
		right_value = rhsvalue(state(k,z,xxi,p),K[right_ind],i_z,i_xxi,right_ind, EV, p, trash);
		if (left_value>right_value) {
			Vplus[index] = left_value;
			koptind[index] = left_ind;
		} else {
			Vplus[index] = right_value;
			koptind[index] = right_ind;
		};
	} else if (right_ind-left_ind==2) {
		double value1 = rhsvalue(state(k,z,xxi,p),K[left_ind],i_z,i_xxi,left_ind, EV, p, trash);
		double value2 = rhsvalue(state(k,z,xxi,p),K[left_ind+1],i_z,i_xxi,left_ind+1, EV, p, trash);
		double value3 = rhsvalue(state(k,z,xxi,p),K[right_ind],i_z,i_xxi,right_ind, EV, p, trash);
		if (value1 < value2) {
			if (value2 < value3) {
				Vplus[index] = value3;
				koptind[index] = right_ind;
			} else {
				Vplus[index] = value2;
				koptind[index] = left_ind+1;
			}
		} else {
			if (value1 < value3) {
				Vplus[index] = value3;
				koptind[index] = right_ind;
			} else { 
				Vplus[index] = value1;
				koptind[index] = left_ind;
			}
		}
	} else {
		int ind1 = left_ind; int ind4 = right_ind;
		int ind2, ind3;
		double value1, value2, value3;
		while (ind4 - ind1 > 2) {
			ind2 = (ind1+ind4)/2;
			ind3 = ind2 + 1;
			value2 = rhsvalue(state(k,z,xxi,p),K[ind2],i_z,i_xxi,ind2, EV, p, trash);
			value3 = rhsvalue(state(k,z,xxi,p),K[ind3],i_z,i_xxi,ind3, EV, p, trash);
			if (value2 < value3) {
				ind1 = ind2;
			} else {
				ind4 = ind3;
			};
		};

		// Now the number of candidates is reduced to three
		value1 = rhsvalue(state(k,z,xxi,p),K[ind1],i_z,i_xxi,ind1, EV, p, trash);
		value2 = rhsvalue(state(k,z,xxi,p),K[ind4-1],i_z,i_xxi,ind4-1, EV, p, trash);
		value3 = rhsvalue(state(k,z,xxi,p),K[ind4],i_z,i_xxi,ind4, EV, p,trash);

		if (value1 < value2) {
			if (value2 < value3) {
				Vplus[index] = value3;
				koptind[index] = ind4;
			} else {
				Vplus[index] = value2;
				koptind[index] = ind4-1;
			}
		} else {
			if (value1 < value3) {
				Vplus[index] = value3;
				koptind[index] = ind4;
			} else { 
				Vplus[index] = value1;
				koptind[index] = ind1;
			}
		}
	}
};

// This functor yields RHS for each (k', k, z). Follwing examples in Thrust doc
struct RHS 
{
	// Data member
	double *K, *Z, *XXI;
	double *V;
	double *Vplus;
	int *koptind;
	double *EV;
	para p;

	// Construct this object, create util from _util, etc.
	__host__ __device__
	RHS(double* K_ptr, double* Z_ptr, double* XXI_ptr,
	double* V_ptr,
	double* Vplus_ptr,
	int* koptind_ptr,
	double* EV_ptr,
	para _p)
	{
		K = K_ptr; Z = Z_ptr; XXI = XXI_ptr;
		V = V_ptr;
		Vplus = Vplus_ptr;
		koptind = koptind_ptr;
		EV = EV_ptr;
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

		// Find and construct state and control, otherwise they won't update in the for loop
		double k =K[i_k]; double z=Z[i_z]; double xxi=XXI[i_xxi];

		// Exploit concavity to update V
		concavemax(k, z, xxi, 0, nk-1, i_k, i_z, i_xxi, K, EV, Vplus, koptind, p);
	};
};	

struct findpolicy {
	// Data member
	double *K, *Z, *XXI,  *EV, *copt, *kopt, *nopt, *mmuopt, *mopt;
	int * koptind;
	para p;

	// Constructor
	__host__ __device__
	findpolicy(double* K_ptr, double* Z_ptr, double* XXI_ptr, int* koptind_ptr, double* EV_ptr, double* copt_ptr, double* kopt_ptr, double* nopt_ptr, double* mmuopt_ptr, double* mopt_ptr, para _p) {
		K = K_ptr;
		Z = Z_ptr;
		XXI = XXI_ptr;
		koptind = koptind_ptr;
		EV = EV_ptr;
		copt = copt_ptr;
		kopt = kopt_ptr;
		nopt = nopt_ptr;
		mmuopt = mmuopt_ptr;
		mopt = mopt_ptr;
		p = _p;
	};

	// Main operator
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

		// Preparation
		double k = K[i_k];
		double z = Z[i_z]; 
		double xxi = XXI[i_xxi];
		double kplus = K[koptind[index]];
		int binding;
		state s(k,z,xxi,p);
		state splus(kplus,z,xxi,p);
		rhsvalue(s, kplus, i_z, i_xxi, koptind[index], EV, p, binding);
		control u;

		// Try not binding first
		if (binding==0) {
			u.compute(s,splus,p,0);
			copt[index] = u.c;
			kopt[index] = splus.k;
			nopt[index] = u.n;
			mmuopt[index] = u.mmu;
			mopt[index] = (1-p.ddelta+(1-u.mmu)*p.ttheta*s.z*pow(s.k,p.ttheta-1)*pow(u.n,1-p.ttheta))/u.c;
		} else {
			u.compute(s,splus,p,1);
			copt[index] = u.c;
			kopt[index] = splus.k;
			nopt[index] = u.n;
			mmuopt[index] = u.mmu;
			mopt[index] = (1-p.ddelta+(1-u.mmu)*p.ttheta*s.z*pow(s.k,p.ttheta-1)*pow(u.n,1-p.ttheta))/u.c;
		};
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
	p.complete(); // complete all implied para, find S-S

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
	host_vector<double> h_V(nk*nz*nxxi, (log(p.css)+p.aalpha*log(1-0.3))/(1-p.bbeta));
	host_vector<double> h_Vplus(nk*nz*nxxi,0);
	host_vector<int> h_koptind(nk*nz*nxxi);
	host_vector<double> h_EV(nk*nz*nxxi,0.0);
	host_vector<double> h_P(nz*nxxi*nz*nxxi, 0);

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
	// guess_linear(h_K, h_Z, h_XXI, h_V1_low, h_V1_high, para, 0.3, 3) ;

	// Copy to the device
	device_vector<double> d_K = h_K;
	device_vector<double> d_Z = h_Z;
	device_vector<double> d_XXI = h_XXI;
	device_vector<double> d_V = h_V;
	device_vector<double> d_Vplus = h_Vplus;
	device_vector<int> d_koptind = h_koptind;
	device_vector<double> d_EV = h_EV;
	device_vector<double> d_P = h_P;

	// Obtain device pointers to be used by cuBLAS
	double* d_K_ptr = raw_pointer_cast(d_K.data());
	double* d_Z_ptr = raw_pointer_cast(d_Z.data());
	double* d_XXI_ptr = raw_pointer_cast(d_XXI.data());
	double* d_V_ptr = raw_pointer_cast(d_V.data());
	double* d_Vplus_ptr = raw_pointer_cast(d_Vplus.data());
	int* d_koptind_ptr = raw_pointer_cast(d_koptind.data());
	double* d_EV_ptr = raw_pointer_cast(d_EV.data());
	double* d_P_ptr = raw_pointer_cast(d_P.data());

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
	
	double diff = 10;  int iter = 0;
	while ((diff>tol)&&(iter<maxiter)){
		// Find EMs for low and high 
		cublasDgemm(handle,
			CUBLAS_OP_N,  
			CUBLAS_OP_T,
			nk, nz*nxxi, nz*nxxi,
			&alpha,
			d_V_ptr, 
			nk, 
			d_P_ptr,
			nz*nxxi,
			&beta,
			d_EV_ptr,
			nk);

		// Directly find the new Value function
		thrust::for_each(
			begin,
			end,
			RHS(d_K_ptr, d_Z_ptr, d_XXI_ptr,
				d_V_ptr,
				d_Vplus_ptr,
				d_koptind_ptr,
				d_EV_ptr,
				p)
		);

		// Find diff 
		diff = transform_reduce(
			make_zip_iterator(make_tuple(d_V.begin(),d_Vplus.begin())),
			make_zip_iterator(make_tuple(d_V.end()  ,d_Vplus.end())),
			myDist(),
			0.0,
			maximum<double>()
			);

		cout << "diff is: "<< diff << endl;
		cout << "Vplus1[117,2,4] (the spike) range is " << d_Vplus[117+2*nk+4*nk*nz] << ", " << d_Vplus[117+2*nk+4*nk*nz] << endl;

		// update correspondence
		d_V = d_Vplus;

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

	// After it converges find the policy
	device_vector<double> d_copt(nk*nz*nxxi);
	device_vector<double> d_kopt(nk*nz*nxxi);
	device_vector<double> d_nopt(nk*nz*nxxi);
	device_vector<double> d_mmuopt(nk*nz*nxxi);
	device_vector<double> d_mopt(nk*nz*nxxi);
	double* d_copt_ptr = raw_pointer_cast(d_copt.data());
	double* d_kopt_ptr = raw_pointer_cast(d_kopt.data());
	double* d_nopt_ptr = raw_pointer_cast(d_nopt.data());
	double* d_mmuopt_ptr = raw_pointer_cast(d_mmuopt.data());
	double* d_mopt_ptr = raw_pointer_cast(d_mopt.data());

	// Find polices
	thrust::for_each(
			begin,
			end,
			findpolicy(d_K_ptr, d_Z_ptr, d_XXI_ptr, d_koptind_ptr, d_EV_ptr, d_copt_ptr, d_kopt_ptr, d_nopt_ptr, d_mmuopt_ptr, d_mopt_ptr,p)
			);

	// Copy back to host and print to file
	h_V = d_V;
	h_EV = d_EV;
	h_koptind = d_koptind;
	
	// Save the decision variables
	host_vector<double> h_copt = d_copt;
	host_vector<double> h_kopt = d_kopt;
	host_vector<double> h_nopt = d_nopt;
	host_vector<double> h_mmuopt = d_mmuopt;
	host_vector<double> h_mopt = d_mopt;

	save_vec(h_K,nk,"./vfi_results/Kgrid.csv");
	save_vec(h_Z,"./vfi_results/Zgrid.csv");
	save_vec(h_XXI,"./vfi_results/XXIgrid.csv");
	save_vec(h_P,"./vfi_results/P.csv");
	save_vec(h_V,"./vfi_results/V.csv");
	save_vec(h_copt,"./vfi_results/copt.csv");
	save_vec(h_kopt,"./vfi_results/kopt.csv");
	save_vec(h_nopt,"./vfi_results/nopt.csv");
	save_vec(h_mmuopt,"./vfi_results/mmuopt.csv");
	save_vec(h_mopt,"./vfi_results/mopt.csv");
	cout << "Policy functions output completed." << endl;

	// Export parameters to MATLAB
	p.exportmatlab("./MATLAB/vfi_para.m");

	return 0;
}
