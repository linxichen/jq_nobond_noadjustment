#define nk 25001
#define nm1 2560
#define nz 13
#define nxxi 13
#define tol 1e-10
#define maxiter 2500
#define kwidth 1.2

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

using namespace std;
using namespace thrust;

// Define an class that contains parameters and steady states
struct para_struct {
	// Model parameters
	double aalpha;
	double bbeta ;
	double ddelta;
	double ttheta;
	double kkappa;
	double ttau  ;
	double xxibar;
	double zbar  ;
	double dbar  ;
	double rrhozz;
	double rrhozxxi;
	double rrhoxxiz;
	double rrhoxxixxi;
	double var_epsz;
	double var_epsxxi;
	double A[4];
	double Ssigma_e[4];

	// Steady States
	double kss;
	double nss;
	double css;
	double wss;
	double dss;
	double mmuss;
	double mkss;
	double yss;


	// Find steady state and find aalpha based steady state target
	__host__ __device__
	void complete() {
		// Fill A and Ssigma_e
		A[0] = rrhozz; A[2] = rrhozxxi;
		A[1] = rrhoxxiz; A[3] = rrhoxxixxi;
		Ssigma_e[0] = var_epsz*var_epsz;
		Ssigma_e[3] = var_epsxxi*var_epsxxi;

		// Find aalpha based on SS computation
		double kovern = pow(xxibar,1/(ttheta-1));
		double covern = pow(kovern,ttheta) - ddelta*kovern;
		mmuss = 1 - ( bbeta*(1-ddelta)-1+xxibar )/( xxibar*(1-bbeta*ttheta)  );
		aalpha = double(0.7/0.3)*(1/covern)*(1-mmuss)*(1-ttheta)*pow(kovern,ttheta);
		double G = ( (1-mmuss)*(1-ttheta)*pow(kovern,ttheta) ) / ( aalpha*covern );
		nss = G/(1+G);
		css = nss*covern;
		kss = nss*kovern;
		wss = aalpha*css/(1-nss);
		dss = css - wss*nss;
		mkss = (1-ddelta+(1-mmuss)*zbar*ttheta*pow(kss,ttheta-1)*pow(nss,1-ttheta))/css;
		yss = zbar*pow(kss,ttheta)*pow(nss,1-ttheta);
	};

	// Export parameters to a .m file in MATLAB syntax
	__host__
	void exportmatlab(std::string filename) {
		std::ofstream fileout(filename.c_str(), std::ofstream::trunc);
		// Accuracy Controls
		fileout << setprecision(16) << "nk=" << nk << ";"<< endl;
		fileout << setprecision(16) << "nz=" << nz << ";"<< endl;
		fileout << setprecision(16) << "nxxi=" << nxxi << ";"<< endl;
		fileout << setprecision(16) << "nm1=" << nm1 << ";"<< endl;

		// Model Parameters
		fileout << setprecision(16) << "aalpha=" << aalpha << ";"<< endl;
		fileout << setprecision(16) << "bbeta=" << bbeta << ";"<< endl;
		fileout << setprecision(16) << "ddelta=" << ddelta << ";"<< endl;
		fileout << setprecision(16) << "ttheta=" << ttheta << ";"<< endl;
		fileout << setprecision(16) << "xxibar=" << xxibar << ";"<< endl;
		fileout << setprecision(16) << "zbar=" << zbar << ";"<< endl;
		fileout << setprecision(16) << "rrhozz=" << rrhozz << ";"<< endl;
		fileout << setprecision(16) << "rrhozxxi=" << rrhozxxi << ";"<< endl;
		fileout << setprecision(16) << "rrhoxxiz=" << rrhoxxiz << ";"<< endl;
		fileout << setprecision(16) << "rrhoxxixxi=" << rrhoxxixxi << ";"<< endl;
		fileout << setprecision(16) << "ssigmaepsz=" << sqrt(var_epsz) << ";"<< endl;
		fileout << setprecision(16) << "ssigmaepsxxi=" << sqrt(var_epsxxi) << ";"<< endl;

		// Steady States
		fileout << setprecision(16) << "kss=" << kss << ";"<< endl;
		fileout << setprecision(16) << "nss=" << nss << ";"<< endl;
		fileout << setprecision(16) << "css=" << css << ";"<< endl;
		fileout << setprecision(16) << "wss=" << wss << ";"<< endl;
		fileout << setprecision(16) << "dss=" << dss << ";"<< endl;
		fileout << setprecision(16) << "mmuss=" << mmuss << ";"<< endl;
		fileout << setprecision(16) << "mkss=" << mkss << ";"<< endl;
		fileout << setprecision(16) << "yss=" << yss << ";"<< endl;
		fileout.close();
	};
};

void guess_linear(const host_vector<double> K, const host_vector<double> Z, const host_vector<double> XXI, host_vector<double> & V1_low, host_vector<double> & V1_high, para_struct para, double factor_low, double factor_high) {
	// Initialize matrices
	int n = 9; int n_jump = 8; int n_shock = 2;
	host_vector<double> A(n*n,0); 
	host_vector<double> B(n*n,0); 
	host_vector<double> C(n*n_shock,0); 
	host_vector<double> rrho(n_shock*n_shock,0);
   	host_vector<double> Pphi(n*(n-n_jump+n_shock),0);

	// Fill in matrices.
	// HH Budget. Correct.
	B[0+3*n] = para.nss;
	B[0+2*n] = para.wss;
	B[0+4*n] = 1;
	B[0+1*n] = -1;

	// Labor Demand. Correct
	B[1+5*n] = (para.ttheta-1)*para.yss/para.nss;
	B[1+6*n] = (1-para.ttheta)*(1-para.mmuss)/para.nss;
	B[1+2*n] = -(1-para.ttheta)*(1-para.mmuss)*para.yss/(para.nss*para.nss);
	B[1+3*n] = -1;

	// Labor Supply. Correct
	B[2+1*n] = para.aalpha/(1-para.nss);
	B[2+2*n] = para.aalpha*para.css/((1-para.nss)*(1-para.nss));
	B[2+3*n] = -1;

	// Capital Demand. Correct.
	A[3+8*n] = para.bbeta; 
	B[3+1*n] = -(1-para.mmuss*para.xxibar)/(para.css*para.css); 
	B[3+5*n] = -para.xxibar/para.css; 
	C[3+1*n] = -para.mmuss*para.xxibar/para.css;

	// Resource Constraint. Correct
	A[4+0*n] = 1; 
	B[4+0*n] = 1-para.ddelta; 
	B[4+6*n] = 1; 
	B[4+1*n] = -1;

	// Financial Constraint. Fixed.
	A[5+0*n] = para.xxibar;
	B[5+6*n] = 1;
	C[5+1*n] = -para.xxibar*para.kss;

	// Output Definition. Correct
	C[6+0*n] = para.yss;
	B[6+0*n] = para.ttheta*para.yss/para.kss;
	B[6+2*n] = (1-para.ttheta)*para.yss/para.nss;
	B[6+6*n] = -1;

	// Investment Definition. Correct
	A[7+0*n] = 1;
	B[7+7*n] = 1;
	B[7+0*n] = 1-para.ddelta;

	// MK defintion:
	B[8+1*n] = -pow(para.css,-2)*(1-para.ddelta+(1-para.mmuss)*para.ttheta*para.yss/para.kss); 
	B[8+5*n] = -para.ttheta*para.yss/(para.css*para.kss); 
	B[8+6*n] = (1-para.mmuss)*para.ttheta/(para.css*para.kss); 
	B[8+0*n] = -(1-para.mmuss)*para.ttheta*para.yss*pow(para.kss,-2)/para.css;
	B[8+8*n] = -1;

	for (int i=0; i< n_shock*n_shock; i++) {
		rrho[i] = para.A[i];
	};

	// Call linear solver
	linearQZ(A.data(),B.data(),C.data(),rrho.data(),n,n_jump,n_shock,Pphi.data());

	// Create guesses.
	for (int i_k=0; i_k<nk; i_k++) {
		for (int i_z = 0; i_z < nz; i_z++) {
			for (int i_xxi = 0; i_xxi < nxxi; i_xxi++) {
				double temp = para.mkss+Pphi[8+0*9]*(K[i_k]-para.kss) + Pphi[8+1*9]*(log(Z[i_z])-log(para.zbar))+ Pphi[8+2*9]*(log(XXI[i_xxi])-log(para.xxibar));
				V1_low[i_k+nk*i_z+nk*nz*i_xxi] = factor_low*temp;
				V1_high[i_k+nk*i_z+nk*nz*i_xxi] = factor_high*temp;
			};
		};
	};
};

// Define state struct that contain things known to shrink functor
struct state_struct {
	// Data member
	double k, z, xxi, zkttheta, kplus;

	// Constructor
	__host__ __device__
	state_struct(double _k, double _z, double _xxi, double _zkttheta, double _kplus) {
		k = _k;
		z = _z;
		xxi = _xxi;
		zkttheta = _zkttheta;
		kplus = _kplus;
	};
	// Pack Things into
	__host__ __device__
	void load(double _k, double _z, double _xxi, double _zkttheta, double _kplus) {
		k = _k;
		z = _z;
		xxi = _xxi;
		zkttheta = _zkttheta;
		kplus = _kplus;
	};

	// Replace only the kplus
	__host__ __device__
	void replace(double _kplus) {
		kplus = _kplus;
	};
};

struct case2_hour {
	// Data Member
	double c0, coneminusttheta, cminusttheta;
	double ttheta;

	// Constructor
	__host__ __device__
	case2_hour(state_struct s, para_struct para) {
		c0 = para.aalpha*(1-para.ddelta)*s.k - para.aalpha*s.kplus;
		coneminusttheta = (para.aalpha+1-para.ttheta)*s.zkttheta;
		cminusttheta = (para.ttheta-1)*s.zkttheta;
		ttheta = para.ttheta;
	};

	// The function
	__host__ __device__
	double operator()(double n) {
		return c0 + coneminusttheta*pow(n,1-ttheta) + cminusttheta*pow(n,-ttheta);
	};

	// The derivative
	__host__ __device__
	double prime(double n) {
		return (1-ttheta)*coneminusttheta*pow(n,-ttheta) + (-ttheta)*cminusttheta*pow(n,-ttheta-1);
	};
};

// Define a struct that contains variables implied by augmented state
struct control_struct {
	// Data member
	double c, n, mmu, Y;

	// Constructor and finding the control variables
	__host__ __device__
	void compute(state_struct state, para_struct para, int binding) {
		if (binding == 1) {
			// Case 1: Binding
			double k = state.k;
			double xxi = state.xxi;
			double zkttheta = state.zkttheta;
			double kplus = state.kplus;
			n = pow(xxi*kplus/(zkttheta),1/(1-para.ttheta));
			Y = xxi*kplus;
			c = Y+(1-para.ddelta)*k-kplus;
			mmu = 1-para.aalpha*c/((1-n)*(1-para.ttheta)*Y/n);	
		};

		if (binding == 0) {
			// Case 2: Not Binding
			double k = state.k;
			double zkttheta = state.zkttheta;
			double kplus = state.kplus;
			n = newton(case2_hour(state,para),1e-5,1.0-1e-5,0.3);
			Y = zkttheta*pow(n,1-para.ttheta);
			mmu = 0;
			c = Y+(1-para.ddelta)*k-kplus;
		};
	};
};

/// This function finds the value of RHS given k', k, z, xxi
__host__ __device__
double rhsvalue (state_struct s, int i_z, int i_xxi, int i_kplus, double* EV, para_struct para, int & binding) {
	double temp1 = -9999999;
	double temp2 = -9999999;
	control_struct u1, u2;
	u1.compute(s,para,1);
	u2.compute(s,para,0);
	if (
			// (u1.mmu >= 0) &&
			(u1.c > 0 ) &&
			(u1.n > 0 ) &&
			(u1.n < 1 )
	   )
	{
		temp1 = log(u1.c) + para.aalpha*log(1-u1.n) + para.bbeta*EV[i_kplus+i_z*nk+i_xxi*nk*nz];
	};
	if (
			(s.xxi*s.kplus > u2.Y) &&
			(u2.c > 0 ) &&
			(u2.n > 0 ) &&
			(u2.n < 1 )
	   )
	{
		temp2 = log(u2.c) + para.aalpha*log(1-u2.n) + para.bbeta*EV[i_kplus+i_z*nk+i_xxi*nk*nz];
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
void concavemax(double k, double z, double xxi, double zkttheta, const int left_ind, const int right_ind, const int i_k,const int i_z, const int i_xxi, 
		double* K, double* EV, double* Vplus, int* koptind, para_struct para) {
	int index = i_k + i_z*nk + i_xxi*nk*nz;
	int trash;

	if (right_ind-left_ind==1) {
		double left_value, right_value;
		left_value = rhsvalue(state_struct(k,z,xxi,zkttheta,K[left_ind]),i_z,i_xxi,left_ind, EV, para, trash);
		right_value = rhsvalue(state_struct(k,z,xxi,zkttheta,K[right_ind]),i_z,i_xxi,right_ind, EV, para, trash);
		if (left_value>right_value) {
			Vplus[index] = left_value;
			koptind[index] = left_ind;
		} else {
			Vplus[index] = right_value;
			koptind[index] = right_ind;
		};
	} else if (right_ind-left_ind==2) {
		double value1 = rhsvalue(state_struct(k,z,xxi,zkttheta,K[left_ind]),i_z,i_xxi,left_ind, EV, para, trash);
		double value2 = rhsvalue(state_struct(k,z,xxi,zkttheta,K[left_ind+1]),i_z,i_xxi,left_ind+1, EV, para, trash);
		double value3 = rhsvalue(state_struct(k,z,xxi,zkttheta,K[right_ind]),i_z,i_xxi,right_ind, EV, para, trash);
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
			value2 = rhsvalue(state_struct(k,z,xxi,zkttheta,K[ind2]),i_z,i_xxi,ind2, EV, para, trash);
			value3 = rhsvalue(state_struct(k,z,xxi,zkttheta,K[ind3]),i_z,i_xxi,ind3, EV, para, trash);
			if (value2 < value3) {
				ind1 = ind2;
			} else {
				ind4 = ind3;
			};
		};

		// Now the number of candidates is reduced to three
		value1 = rhsvalue(state_struct(k,z,xxi,zkttheta,K[ind1]),i_z,i_xxi,ind1, EV, para, trash);
		value2 = rhsvalue(state_struct(k,z,xxi,zkttheta,K[ind4-1]),i_z,i_xxi,ind4-1, EV, para, trash);
		value3 = rhsvalue(state_struct(k,z,xxi,zkttheta,K[ind4]),i_z,i_xxi,ind4, EV, para,trash);

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
	para_struct para;

	// Construct this object, create util from _util, etc.
	__host__ __device__
	RHS(double* K_ptr, double* Z_ptr, double* XXI_ptr,
	double* V_ptr,
	double* Vplus_ptr,
	int* koptind_ptr,
	double* EV_ptr,
	para_struct _para)
	{
		K = K_ptr; Z = Z_ptr; XXI = XXI_ptr;
		V = V_ptr;
		Vplus = Vplus_ptr;
		koptind = koptind_ptr;
		EV = EV_ptr;
		para = _para;
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
		double zkttheta = z*pow(k,para.ttheta); 

		// Exploit concavity to update V
		concavemax(k, z, xxi, zkttheta, 0, nk-1, i_k, i_z, i_xxi, K, EV, Vplus, koptind, para);
	};
};	

struct findpolicy {
	// Data member
	double *K, *Z, *XXI,  *EV, *copt, *kopt, *nopt, *mmuopt;
	int * koptind;
	para_struct para;

	// Constructor
	__host__ __device__
	findpolicy(double* K_ptr, double* Z_ptr, double* XXI_ptr, int* koptind_ptr, double* EV_ptr, double* copt_ptr, double* kopt_ptr, double* nopt_ptr, double* mmuopt_ptr, para_struct _para) {
		K = K_ptr;
		Z = Z_ptr;
		XXI = XXI_ptr;
		koptind = koptind_ptr;
		EV = EV_ptr;
		copt = copt_ptr;
		kopt = kopt_ptr;
		nopt = nopt_ptr;
		mmuopt = mmuopt_ptr;
		para = _para;
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
		double zkttheta = z*pow(k,para.ttheta);
		int binding;
		state_struct s(k,z,xxi,zkttheta,kplus);
		rhsvalue(s, i_z, i_xxi, koptind[index], EV, para, binding);
		control_struct u;

		// Try not binding first
		if (binding==0) {
			u.compute(s,para,0);
			copt[index] = u.c;
			kopt[index] = s.kplus;
			nopt[index] = u.n;
			mmuopt[index] = u.mmu;
		} else {
			u.compute(s,para,1);
			copt[index] = u.c;
			kopt[index] = s.kplus;
			nopt[index] = u.n;
			mmuopt[index] = u.mmu;
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
	para_struct para;

	// Set Model Parameters
	para.bbeta = 0.9825;
	para.ddelta = 0.025;
	para.ttheta = 0.36;
	para.kkappa = 0.1460;
	para.ttau = 0.3500;
	para.xxibar = 0.11;
	para.zbar = 1.0;
	para.rrhozz = 0.9457;
	para.rrhoxxiz = 0.0321;
	para.rrhozxxi =-0.0091;
	para.rrhoxxixxi = 0.9703;
	para.var_epsz = 0.0045*0.0045;
	para.var_epsxxi = 0.0098*0.0098;
	para.complete(); // complete all implied para, find S-S

	cout << setprecision(16) << "kss: " << para.kss << endl;
	cout << setprecision(16) << "zss: " << para.zbar << endl;
	cout << setprecision(16) << "xxiss: " <<para.xxibar << endl;
	cout << setprecision(16) << "mkss: " << para.mkss << endl;
	cout << setprecision(16) << "dss: " << para.dss << endl;
	cout << setprecision(16) << "css: " << para.css << endl;
	cout << setprecision(16) << "nss: " << para.nss << endl;
	cout << setprecision(16) << "wss: " << para.wss << endl;
	cout << setprecision(16) << "mmuss: " << para.mmuss << endl;
	cout << setprecision(16) << "aalpha: " << para.aalpha << endl;
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
	host_vector<double> h_V(nk*nz*nxxi, (log(para.css)+para.aalpha*log(1-0.3))/(1-para.bbeta));
	host_vector<double> h_Vplus(nk*nz*nxxi,0);
	host_vector<int> h_koptind(nk*nz*nxxi);
	host_vector<double> h_EV(nk*nz*nxxi,0.0);
	host_vector<double> h_P(nz*nxxi*nz*nxxi, 0);

	// Create capital grid
	double minK = 1/kwidth*para.kss;
	double maxK = kwidth*para.kss;
	linspace(minK,maxK,nk,raw_pointer_cast(h_K.data()));

	// Create shocks grids
	host_vector<double> h_shockgrids(2*nz);
	double* h_shockgrids_ptr = raw_pointer_cast(h_shockgrids.data());
	double* h_P_ptr = raw_pointer_cast(h_P.data());
	gridgen_fptr linspace_fptr = &linspace; // select linspace as grid gen
	tauchen_vec(2,nz,3,para.A,para.Ssigma_e,h_shockgrids_ptr,h_P_ptr,linspace_fptr);
	for (int i_shock = 0; i_shock < nz; i_shock++) {
		h_Z[i_shock] = para.zbar*exp(h_shockgrids[i_shock+0*nz]);
		h_XXI[i_shock] = para.xxibar*exp(h_shockgrids[i_shock+1*nz]);
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
				para)
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
	double* d_copt_ptr = raw_pointer_cast(d_copt.data());
	double* d_kopt_ptr = raw_pointer_cast(d_kopt.data());
	double* d_nopt_ptr = raw_pointer_cast(d_nopt.data());
	double* d_mmuopt_ptr = raw_pointer_cast(d_mmuopt.data());

	// Find polices
	thrust::for_each(
			begin,
			end,
			findpolicy(d_K_ptr, d_Z_ptr, d_XXI_ptr, d_koptind_ptr, d_EV_ptr, d_copt_ptr, d_kopt_ptr, d_nopt_ptr, d_mmuopt_ptr, para)
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
	// host_vector<double> h_dopt(nk*nz*nxxi);
	// host_vector<double> h_wopt(nk*nz*nxxi);

	save_vec(h_K,nk,"./vfi_results/Kgrid.csv");
	save_vec(h_Z,"./vfi_results/Zgrid.csv");
	save_vec(h_XXI,"./vfi_results/XXIgrid.csv");
	save_vec(h_P,"./vfi_results/P.csv");
	save_vec(h_V,"./vfi_results/V.csv");
	save_vec(h_copt,"./vfi_results/copt.csv");
	save_vec(h_kopt,"./vfi_results/kopt.csv");
	save_vec(h_nopt,"./vfi_results/nopt.csv");
	save_vec(h_mmuopt,"./vfi_results/mmuopt.csv");
	cout << "Policy functions output completed." << endl;

	// Export parameters to MATLAB
	para.exportmatlab("./MATLAB/vfi_para.m");

	// Compute the Euler equation error
	double eee = -999999;
	for (int i_k = 0; i_k < nk; i_k++) {
		for (int i_z = 0; i_z < nz; i_z++) {
			for (int i_xxi = 0; i_xxi < nxxi; i_xxi++) {
				int index = i_k + i_z*nk + i_xxi*nk*nz;
				double ctoday = h_copt[index];
				int i_kplus = h_koptind[index];
				double sum = 0;
				for (int i_zplus = 0; i_zplus < nz; i_zplus++) {
					for (int i_xxiplus = 0; i_xxiplus < nxxi; i_xxiplus++) {
						int indexplus = i_kplus+i_zplus*nk+i_xxiplus*nk*nz;
						sum += para.bbeta*h_P[i_z + i_xxi*nz + i_zplus*nz*nxxi + i_xxiplus*nz*nxxi*nz]*(1 - para.ddelta + (1-h_mmuopt[indexplus])*para.ttheta*h_Z[i_zplus]*pow(h_K[i_kplus],para.ttheta-1)*pow(h_nopt[indexplus],1-para.ttheta))/h_copt[indexplus];
					};
				};
				double ctmr = (1-h_mmuopt[index]*h_XXI[i_xxi])/sum;
				double eee_temp = abs(ctoday/ctmr-1);
				eee = max(eee, eee_temp);
			};
		};
	};

	cout << "Euler equation error = " << eee << endl;
	return 0;
}
