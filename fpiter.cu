/* Includes, system */
#include <iostream>
#include <iomanip>
#include <fstream>

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
#include <cublas_v2.h>
#include "cuda_helpers.h"

// Includes, Armadillo
#include "cppcode.h"

// Includes, model specific
#include "model.h"

using namespace std;
using namespace thrust;

// #define M_PI 3.14159265358979323846264338328
#define nk 12 
#define nb 1 
#define nz 13 
#define nxxi 13 
#define nm1 501 
#define pk 6
#define pz 6
#define pxxi 6
#define tol 1e-7
#define maxiter 1000
#define kwidth 1.2
#define bwidth 1.15 
#define llambda 0.95

struct findpolicy {
	// Data member
	double *K, *Z, *XXI,  *coeff, *copt, *kopt, *nopt, *mmuopt ;
	double *K_cheby, *Z_cheby, *XXI_cheby;
	int nkout;
	para p;

	// Constructor
	__host__ __device__
	findpolicy(double* K_ptr, double* K_cheby_ptr, double* Z_ptr, double* Z_cheby_ptr, double* XXI_ptr, double* XXI_cheby_ptr, double* coeff_ptr, double* copt_ptr, double* kopt_ptr, double* nopt_ptr, double* mmuopt_ptr, int _nkout, para _p) {
		K = K_ptr;
		Z = Z_ptr;
		XXI = XXI_ptr;
		K_cheby = K_cheby_ptr;
		Z_cheby = Z_cheby_ptr;
		XXI_cheby = XXI_cheby_ptr;
		coeff = coeff_ptr;
		copt = copt_ptr;
		kopt = kopt_ptr;
		nopt = nopt_ptr;
		mmuopt = mmuopt_ptr;
		nkout = _nkout;
		p = _p;
	};

	// Main operator
	__host__ __device__
	void operator()(int index) {
		// Perform ind2sub
		int subs[3];
		int size_vec[3];
		size_vec[0] = nkout;
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
		state s(k,z,xxi,p);

		// Find the current M
		double k_cheby = K_cheby[i_k];
		double z_cheby = Z_cheby[i_z];
		double xxi_cheby = XXI_cheby[i_xxi];
		double arg[3]; 
		arg[0] = k_cheby;
		arg[1] = z_cheby;
		arg[2] = xxi_cheby;
		size_vec[0] = pk+1;
		size_vec[1] = pz+1;
		size_vec[2] = pxxi+1;
		int temp_subs[3];
		double m1 = chebyeval_multi(3,arg,size_vec,temp_subs,coeff);

		// Try not binding first
		control u;
		shadow m(m1);
		u.compute(s,m,p,0);
		if (s.xxi*u.kplus>u.Y) {
			copt[index] = u.c;
			kopt[index] = u.kplus;
			nopt[index] = u.n;
			mmuopt[index] = u.mmu;
		} else {
			u.compute(s,m,p,1);
			copt[index] = u.c;
			kopt[index] = u.kplus;
			nopt[index] = u.n;
			mmuopt[index] = u.mmu;
		};
	};
};

void guess_vfi(const host_vector<double> K, const host_vector<double> Z, const host_vector<double> XXI, host_vector<double> & M, para p, double factor) {
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
				M[i_k+nk*i_z+nk*nz*i_xxi] = factor*temp;
			};
		};
	};
};

void guess_linear(const host_vector<double> K, const host_vector<double> Z, const host_vector<double> XXI, host_vector<double> & M, para p, double factor) {
	// Initialize matrices
	int n = 9; int n_jump = 8; int n_shock = 2;
	host_vector<double> A(n*n,0); 
	host_vector<double> B(n*n,0); 
	host_vector<double> C(n*n_shock,0); 
	host_vector<double> rrho(n_shock*n_shock,0);
   	host_vector<double> Pphi(n*(n-n_jump+n_shock),0);

	// Fill in matrices.
	linearizedmodel(A.data(), B.data(), C.data(), rrho.data(), n, n_shock, p);

	// Call linear solver
	linearQZ(A.data(),B.data(),C.data(),rrho.data(),n,n_jump,n_shock,Pphi.data());

	// Create guesses.
	for (int i_k=0; i_k<nk; i_k++) {
		for (int i_z = 0; i_z < nz; i_z++) {
			for (int i_xxi = 0; i_xxi < nxxi; i_xxi++) {
				double temp = p.mkss+Pphi[8+0*9]*(K[i_k]-p.kss) + Pphi[8+1*9]*(log(Z[i_z])-log(p.zbar))+ Pphi[8+2*9]*(log(XXI[i_xxi])-log(p.xxibar));
				M[i_k+nk*i_z+nk*nz*i_xxi] = factor*temp;
			};
		};
	};
};

// This functor find new M at each state
struct findnewM
{
	// Data Member
	double *K, *K_cheby;
	double *Z, *Z_cheby;
	double *XXI, *XXI_cheby;
	double *P;
	double *M;
	double *M_new;
	double *coeff;
	double minK, maxK;
	para p;

	// Construct this object, create util from _util, etc.
	findnewM(double* K_ptr, double* K_cheby_ptr, double* Z_ptr, double* Z_cheby_ptr, double*XXI_ptr, double* XXI_cheby_ptr,
			 double* P_ptr, double* M_ptr, double* M_new_ptr, double* coeff_ptr, double _minK, double _maxK, para _p)
	{
		K = K_ptr; K_cheby = K_cheby_ptr;
		Z = Z_ptr; Z_cheby = Z_cheby_ptr;
		XXI = XXI_ptr; XXI_cheby = XXI_cheby_ptr;
		P = P_ptr;
		M = M_ptr; M_new = M_new_ptr;
		coeff = coeff_ptr;
		minK = _minK; maxK = _maxK;
		p = _p;
	};

	__host__ __device__
	void operator()(int index) {
		int i_xxi = index/(nk*nz);
		int i_z = (index - i_xxi*nk*nz)/(nk);
		int i_k = index - i_xxi*nk*nz - i_z*nk;
		double k, z, xxi;
		double kplus_cheby, zplus_cheby, xxiplus_cheby;	
		double ctilde;
		double EM;
		double ttheta = p.ttheta;
		double ddelta = p.ddelta;
		double bbeta = p.bbeta;

		// Load Variables
		k = K[i_k]; z = Z[i_z]; xxi = XXI[i_xxi];
		state s(K[i_k],Z[i_z],XXI[i_xxi],p); shadow m(M[index]);
		control u1, u2;

		// Case 1: Binding
		u1.compute(s,m,p,1);
		kplus_cheby = -1 + (u1.kplus-minK)/(maxK-minK)*(2);
		// printf("kplus_cheby=%f\n",kplus_cheby);

		EM = 0;
		for (int i_zplus=0; i_zplus<nz; i_zplus++) {
			zplus_cheby = Z_cheby[i_zplus];
			for (int i_xxiplus=0; i_xxiplus<nxxi; i_xxiplus++) {
				xxiplus_cheby = XXI_cheby[i_xxiplus];
				double arg[3]; 
				arg[0] = kplus_cheby;
				arg[1] = zplus_cheby;
				arg[2] = xxiplus_cheby;
				int size_vec[3];
				size_vec[0] = pk+1;
				size_vec[1] = pz+1;
				size_vec[2] = pxxi+1;
				int temp_subs[3];
				EM += P[i_z+i_xxi*nz+nz*nxxi*i_zplus+nz*nxxi*nz*i_xxiplus]*chebyeval_multi(3,arg,size_vec,temp_subs,coeff);
				// EM += P[i_z+i_xxi*nz+nz*nxxi*i_zplus+nz*nxxi*nz*i_xxiplus]*chebyeval_multi_old(kplus_cheby,zplus_cheby,xxiplus_cheby,coeff);
			};
		};
		ctilde = (1-u1.mmu*s.xxi)/(bbeta*EM);
		// printf("case1 EM = %f\n",EM);

		// Check whether implied policy functions make sense
		if (
			(u1.c>0) && (ctilde>0) && (u1.mmu>=0) && (u1.kplus > minK) && (u1.kplus<maxK) && (u1.n>0) && (u1.n<1)
		   )
		{
			M_new[index] = ((1-u1.mmu)*z*ttheta*pow(k,ttheta-1)*pow(u1.n,1-ttheta)+1-ddelta)/ctilde;
			return;
		} else {
			goto case2;
		}; 

		// Case 2: Not Binding
		case2:
		u2.compute(s,m,p,0);
		kplus_cheby = -1 + (u2.kplus-minK)/(maxK-minK)*(2);

		EM = 0;
		for (int i_zplus=0; i_zplus<nz; i_zplus++) {
			zplus_cheby = Z_cheby[i_zplus];
			for (int i_xxiplus=0; i_xxiplus<nxxi; i_xxiplus++) {
				xxiplus_cheby = XXI_cheby[i_xxiplus];
				double arg[3]; 
				arg[0] = kplus_cheby;
				arg[1] = zplus_cheby;
				arg[2] = xxiplus_cheby;
				int size_vec[3];
				size_vec[0] = pk+1;
				size_vec[1] = pz+1;
				size_vec[2] = pxxi+1;
				int temp_subs[3];
				EM += P[i_z+i_xxi*nz+nz*nxxi*i_zplus+nz*nxxi*nz*i_xxiplus]*chebyeval_multi(3,arg,size_vec,temp_subs,coeff);
				// EM += P[i_z+i_xxi*nz+nz*nxxi*i_zplus+nz*nxxi*nz*i_xxiplus]*chebyeval_multi_old(kplus_cheby,zplus_cheby,xxiplus_cheby,coeff);
			};
		};
		ctilde = (1-u2.mmu*s.xxi)/(bbeta*EM);
		// printf("case2 EM = %f\n",EM);

		if (
			(ctilde>0) && (u2.c>0) && (u2.kplus>minK) && (u2.kplus<maxK) && (u2.kplus*xxi>u2.Y) && (u2.n>0) && (u2.n<1)
		   )
		{
			M_new[index] = ((1-u2.mmu)*z*ttheta*pow(k,ttheta-1)*pow(u2.n,1-ttheta)+1-ddelta)/ctilde;
		} else {
			printf("No solution at k=%f, z=%f, xxi=%f, m=%f.\nPolicies are: c=%f, ctilde=%f, kplus=%f, mu=%f, n=%f\n====================================================\n",s.k,s.z,s.xxi,m.m1,u2.c,ctilde,u2.kplus,u2.mmu,u2.n);
		};
	};
};

// This functor yields X, storing basis function values at each state. See Der Hann notes
struct findbasis
{
	// Data member
	double *K_cheby, *Z_cheby, *XXI_cheby;
	double *X_cheby;
	// double *V1_low;
	// double *V1_high;
	// double *Vplus1_low;
	// double *Vplus1_high;
	// double *EM1_low;
	// double *EM1_high;
	// double *flag;

	// Construct this object, create util from _util, etc.
	__host__ __device__
	findbasis(double* K_cheby_ptr, double* Z_cheby_ptr, double* XXI_cheby_ptr, double* X_cheby_ptr)
	// double* V1_low_ptr,
	// double* V1_high_ptr,
	// double* Vplus1_low_ptr,
	// double* Vplus1_high_ptr,
	// double* EM1_low_ptr,
	// double* EM1_high_ptr,
	// double* flag_ptr)
	{
		K_cheby = K_cheby_ptr; Z_cheby = Z_cheby_ptr; XXI_cheby = XXI_cheby_ptr;
		X_cheby = X_cheby_ptr;
		// V1_low = V1_low_ptr;
		// V1_high = V1_high_ptr;
		// Vplus1_low = Vplus1_low_ptr;
		// Vplus1_high = Vplus1_high_ptr;
		// EM1_low = EM1_low_ptr;
		// EM1_high = EM1_high_ptr;
		// flag = flag_ptr;
	};

	__host__ __device__
	void operator()(int index) {
		int j_xxi =index/(nk*nz*nxxi*(1+pk)*(1+pz)); // j stands for the order of cheby polynomial
		int j_z = (index - j_xxi*nk*nz*nxxi*(1+pk)*(1+pz)) / (nk*nz*nxxi*(1+pk));
		int j_k = (index - j_xxi*nk*nz*nxxi*(1+pk)*(1+pz) - j_z*nk*nz*nxxi*(1+pk)) / (nk*nz*nxxi);
		int i_xxi=(index - j_xxi*nk*nz*nxxi*(1+pk)*(1+pz) - j_z*nk*nz*nxxi*(1+pk) - j_k*nk*nz*nxxi) / (nk*nz);
		int i_z = (index - j_xxi*nk*nz*nxxi*(1+pk)*(1+pz) - j_z*nk*nz*nxxi*(1+pk) - j_k*nk*nz*nxxi - i_xxi*nk*nz) / (nk);
		int i_k = (index - j_xxi*nk*nz*nxxi*(1+pk)*(1+pz) - j_z*nk*nz*nxxi*(1+pk) - j_k*nk*nz*nxxi - i_xxi*nk*nz - i_z*nk) / (1);

		double k =K_cheby[i_k]; double z=Z_cheby[i_z]; double xxi=XXI_cheby[i_xxi];
		X_cheby[index] = chebypoly(j_k,k)*chebypoly(j_z,z)*chebypoly(j_xxi,xxi);
		// printf("Coordinates=(%i, %i, %i, %i, %i, %i), Value=%f\n", i_k, i_z, i_xxi, j_k, j_z, j_xxi, X_cheby[index]);
	};
};	

// This functor find the dampened coefficient vector
struct updatecoeff 
{
	// Data Member
	double *coeff, *coeff_temp, *coeff_new;

	// Constructor
	updatecoeff(double* coeff_ptr, double* coeff_temp_ptr, double* coeff_new_ptr ) {
		coeff = coeff_ptr; coeff_temp = coeff_temp_ptr; coeff_new = coeff_new_ptr;
	};

	__host__ __device__
	void operator()(int index) {
		coeff_new[index] = llambda*coeff_temp[index] + (1-llambda)*coeff[index];
	};
};

// This functor calculates the error
struct myMinus {
	// Tuple is (V1low,Vplus1low,V1high,Vplus1high,...)
	template <typename Tuple>
	__host__ __device__
	double operator()(Tuple t)
	{
		return abs(get<0>(t)-get<1>(t)) ;
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

// Main
int main(int argc, char** argv)
{
	// Testing Ground
	double x[2]; x[0] = 0.3; x[1] = 0.4;
	int size_vec[2]; size_vec[0] = 3; size_vec[1] = 3;
	int subs[2];
	double coeff[9];
	coeff[0] = 1;
	coeff[1] = 2;
	coeff[2] = 3;
	coeff[3] = 4;
	coeff[4] = 5;
	coeff[5] = 6;
	coeff[6] = 7;
	coeff[7] = 8;
	coeff[8] = 9;
	cout << chebyeval_multi(2,x,size_vec,subs,coeff) << endl;

	// Set Model Parameters
	para p;
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
	host_vector<double> h_K_cheby(nk); 
	host_vector<double> h_Z(nz);
	host_vector<double> h_Z_cheby(nz);
	host_vector<double> h_XXI(nxxi);
	host_vector<double> h_XXI_cheby(nxxi);
	host_vector<double> h_P(nz*nxxi*nz*nxxi);
	host_vector<double> h_flag(nk*nb*nz*nxxi, 0); 
	host_vector<double> h_X(nk*nz*nxxi*(1+pk)*(1+pz)*(1+pxxi)); 
	host_vector<double> h_projector((1+pk)*(1+pz)*(1+pxxi)*nk*nz*nxxi); 
	host_vector<double> h_M(nk*nz*nxxi,0.5*p.mkss); 
	host_vector<double> h_M_new(nk*nz*nxxi,0.5*p.mkss); 
	host_vector<double> h_coeff((1+pk)*(1+pz)*(1+pxxi),0.1); 
	host_vector<double> h_coeff_temp((1+pk)*(1+pz)*(1+pxxi),0.1); 
	host_vector<double> h_coeff_new((1+pk)*(1+pz)*(1+pxxi),0.1); 
	
	// Create capital grid
	double* h_K_cheby_ptr = raw_pointer_cast(h_K_cheby.data());
	chebyroots(nk,h_K_cheby_ptr);
	h_K = h_K_cheby;
	double* h_K_ptr = raw_pointer_cast(h_K.data());
	double minK = (1/kwidth)*p.kss;
	double maxK = (1*kwidth)*p.kss;
	cout << "minK: " << minK << endl;
	cout << "maxK: " << maxK << endl;
	fromchebydomain(minK, maxK, nk, h_K_ptr);

	// Create shocks grids
	host_vector<double> h_shockgrids(2*nz);
	double* h_shockgrids_ptr = raw_pointer_cast(h_shockgrids.data());
	double* h_P_ptr = raw_pointer_cast(h_P.data());
	gridgen_fptr chebyspace_fptr = &chebyspace; // select linspace as grid gen
	tauchen_vec(2,nz,4,p.A,p.Ssigma_e,h_shockgrids_ptr,h_P_ptr,chebyspace_fptr);
	for (int i_shock = 0; i_shock < nz; i_shock++) {
		h_Z[i_shock] = p.zbar*exp(h_shockgrids[i_shock+0*nz]);
		h_XXI[i_shock] = p.xxibar*exp(h_shockgrids[i_shock+1*nz]);
	};
	double* h_Z_cheby_ptr = raw_pointer_cast(h_Z_cheby.data());
	double* h_XXI_cheby_ptr = raw_pointer_cast(h_XXI_cheby.data());
	chebyroots(nz,h_Z_cheby_ptr);
	chebyroots(nxxi,h_XXI_cheby_ptr);

	// Create Initial M generated from linear solution
	guess_linear(h_K,h_Z, h_XXI, h_M, p, 1.0);
	// guess_vfi(h_K, h_Z, h_XXI, h_M, p, 1.0); 
	save_vec(h_M,"./fpiter_results/M_guess.csv");

	// Copy to the device
	device_vector<double> d_K = h_K;
	device_vector<double> d_Z = h_Z;
	device_vector<double> d_XXI = h_XXI;

	device_vector<double> d_K_cheby = h_K_cheby;
	device_vector<double> d_Z_cheby = h_Z_cheby;
	device_vector<double> d_XXI_cheby = h_XXI_cheby;

	device_vector<double> d_P = h_P;
	device_vector<double> d_flag = h_flag;

	device_vector<double> d_X = h_X;
	device_vector<double> d_projector = h_projector;
	device_vector<double> d_M = h_M;
	device_vector<double> d_M_new = h_M_new;
	device_vector<double> d_coeff = h_coeff;
	device_vector<double> d_coeff_temp = h_coeff_temp;
	device_vector<double> d_coeff_new = h_coeff_new;

	// Obtain device pointers
	double* d_K_ptr = raw_pointer_cast(d_K.data());
	double* d_Z_ptr = raw_pointer_cast(d_Z.data());
	double* d_XXI_ptr = raw_pointer_cast(d_XXI.data());

	double* d_K_cheby_ptr = raw_pointer_cast(d_K_cheby.data());
	double* d_Z_cheby_ptr = raw_pointer_cast(d_Z_cheby.data());
	double* d_XXI_cheby_ptr = raw_pointer_cast(d_XXI_cheby.data());

	double* d_P_ptr = raw_pointer_cast(d_P.data());
	double* d_flag_ptr = raw_pointer_cast(d_flag.data());
	double* d_X_ptr = raw_pointer_cast(d_X.data());

	double* d_M_ptr = raw_pointer_cast(d_M.data());
	double* d_M_new_ptr = raw_pointer_cast(d_M_new.data());
	double* d_projector_ptr = raw_pointer_cast(d_projector.data());
	double* d_coeff_ptr= raw_pointer_cast(d_coeff.data());
	double* d_coeff_temp_ptr= raw_pointer_cast(d_coeff_temp.data());
	double* d_coeff_new_ptr= raw_pointer_cast(d_coeff_new.data());

	// Firstly a virtual index array from 0 to nk*nk*nz
	counting_iterator<int> begin(0);
	counting_iterator<int> end(nk*nz*nxxi*(1+pk)*(1+pz)*(1+pxxi));

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
	
	// Find the projector matrix once and for all
	double* h_X_ptr = raw_pointer_cast(h_X.data());
	double* h_projector_ptr = raw_pointer_cast(h_projector.data());
	thrust::for_each(
		begin,
		end,
		findbasis(d_K_cheby_ptr, d_Z_cheby_ptr, d_XXI_cheby_ptr, d_X_ptr)
	);
	h_X = d_X;
	findprojector(h_X_ptr, nk*nz*nxxi, (1+pk)*(1+pz)*(1+pxxi), h_projector_ptr);
	d_projector = h_projector; // Copy host projector to device

	// Regress M_guess on basis to find initial coefficient
	// This is doing coeff = (X'X)*X'*M
	cublasDgemv(
			handle,	// cuBlas handle
			CUBLAS_OP_N, // N means don't transpose A
			(1+pk)*(1+pz)*(1+pxxi), // # of row in matrix A
			nk*nz*nxxi, // # of col in matrix A
			&alpha, // just 1
			d_projector_ptr, // pointer to matrix A stored in column-major format
			(1+pk)*(1+pz)*(1+pxxi), // leading dimesnion of array to store A, usually # of rows
			d_M_ptr, // pointer to x
			1, // stride of x, usually 1
			&beta, // usually zero
			d_coeff_ptr, // pointer to y
			1 // stride of y
			);

	// Main iterations
	double diff = 10; int iter = 0;
	while ((diff>tol)&&(iter<maxiter)){
		// Find the current M at each state. Does y = A*x/ M = X*coeff
		cublasDgemv(
				handle,	// cuBlas handle
				CUBLAS_OP_N, // N means don't transpose A
				nk*nz*nxxi, // # of row in matrix A
				(1+pk)*(1+pz)*(1+pxxi), // # of col in matrix A
				&alpha, // just 1
				d_X_ptr, // pointer to matrix A stored in column-major format
				nk*nz*nxxi, // leading dimesnion of array to store A, usually # of rows
				d_coeff_ptr, // pointer to x
				1, // stride of x, usually 1
				&beta, // usually zero
				d_M_ptr, // pointer to y
				1 // stride of y
				);

		// Based on current M(k,z,xxi), find implied new M
		thrust::for_each(
				make_counting_iterator(0),
				make_counting_iterator(nk*nz*nxxi),
				findnewM(d_K_ptr, d_K_cheby_ptr, d_Z_ptr, d_Z_cheby_ptr, d_XXI_ptr, d_XXI_cheby_ptr, d_P_ptr, d_M_ptr, d_M_new_ptr, d_coeff_ptr, minK, maxK, p)
				);

		// Regress new M on basis to find temporary coefficient
		cublasDgemv(
				handle,	// cuBlas handle
				CUBLAS_OP_N, // N means don't transpose A
				(1+pk)*(1+pz)*(1+pxxi), // # of row in matrix A
				nk*nz*nxxi, // # of col in matrix A
				&alpha, // just 1
				d_projector_ptr, // pointer to matrix A stored in column-major format
				(1+pk)*(1+pz)*(1+pxxi), // leading dimesnion of array to store A, usually # of rows
				d_M_new_ptr, // pointer to x
				1, // stride of x, usually 1
				&beta, // usually zero
				d_coeff_temp_ptr, // pointer to y
				1 // stride of y
				);

		// Update coefficient with dampening
		thrust::for_each(
			make_counting_iterator(0),
			make_counting_iterator((1+pk)*(1+pz)*(1+pxxi)),
			updatecoeff(d_coeff_ptr,d_coeff_temp_ptr,d_coeff_new_ptr)
		);

		// Compute difference between coefficient vectors
		diff = transform_reduce(
			make_zip_iterator(make_tuple(d_coeff.begin(),d_coeff_new.begin())),
			make_zip_iterator(make_tuple(d_coeff.end(),d_coeff_new.end())),
			myMinus(),
			0.0,
			maximum<double>()
			);

		// Replace old coefficient with new
		d_coeff = d_coeff_new;

		iter++;
		printf("=======================================================\n=== Iteration No. %i finished \n=======================================================\n",iter);
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

	// Find policy
	int nkout = 50001;
	host_vector<double> h_Kgrid(nkout);
	host_vector<double> h_Kgrid_cheby(nkout);
	chebyroots(nkout,h_Kgrid_cheby.data());
	chebyspace(minK,maxK,nkout,h_Kgrid.data());
	device_vector<double> d_Kgrid = h_Kgrid;
	device_vector<double> d_Kgrid_cheby = h_Kgrid_cheby;
	device_vector<double> d_copt(nkout*nz*nxxi);
	device_vector<double> d_kopt(nkout*nz*nxxi);
	device_vector<double> d_nopt(nkout*nz*nxxi);
	device_vector<double> d_mmuopt(nkout*nz*nxxi);
	double* d_Kgrid_ptr = raw_pointer_cast(d_Kgrid.data());
	double* d_Kgrid_cheby_ptr = raw_pointer_cast(d_Kgrid_cheby.data());
	double* d_copt_ptr = raw_pointer_cast(d_copt.data());
	double* d_kopt_ptr = raw_pointer_cast(d_kopt.data());
	double* d_nopt_ptr = raw_pointer_cast(d_nopt.data());
	double* d_mmuopt_ptr = raw_pointer_cast(d_mmuopt.data());

	thrust::for_each(
			make_counting_iterator(0),
			make_counting_iterator(nkout*nz*nxxi),
			findpolicy(d_Kgrid_ptr, d_Kgrid_cheby_ptr, d_Z_ptr, d_Z_cheby_ptr, d_XXI_ptr, d_XXI_cheby_ptr,d_coeff_ptr,d_copt_ptr,d_kopt_ptr,d_nopt_ptr,d_mmuopt_ptr,nkout,p)
			);

	// Copy back to host and print to file
	h_coeff = d_coeff;
	h_M_new = d_M_new;
	h_M = d_M;
	host_vector<double> h_copt = d_copt;
	host_vector<double> h_kopt = d_kopt;
	host_vector<double> h_nopt = d_nopt;
	host_vector<double> h_mmuopt = d_mmuopt;
    save_vec(h_Kgrid,"./fpiter_results/Kgrid.csv");
    save_vec(h_Z,"./fpiter_results/Zgrid.csv");
    save_vec(h_XXI,"./fpiter_results/XXIgrid.csv");
    save_vec(h_P,"./fpiter_results/P.csv");
    save_vec(h_copt,"./fpiter_results/copt.csv");
    save_vec(h_kopt,"./fpiter_results/kopt.csv");
    save_vec(h_nopt,"./fpiter_results/nopt.csv");
    save_vec(h_mmuopt,"./fpiter_results/mmuopt.csv");
    save_vec(h_coeff,"./fpiter_results/coeff.csv");
    save_vec(h_M,"./fpiter_results/M.csv");
    save_vec(h_M_new,"./fpiter_results/M_new.csv");
	p.exportmatlab("./MATLAB/fpiter_para.m");

	// Save accuracy controls
	std::ofstream fileout("./fpiter_results/accuracy.m", std::ofstream::trunc);
	fileout << std::setprecision(16) << "pk=" << pk << ";"<< std::endl;
	fileout << std::setprecision(16) << "pz=" << pz<< ";"<< std::endl;
	fileout << std::setprecision(16) << "pxxi=" << pxxi<< ";"<< std::endl;
	fileout.close();

    return 0;
};


