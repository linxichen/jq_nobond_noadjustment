/* This is a CUDA implementation of Jermann Quadrini 2013 AER  
 * We simulate a equilibrium that collateral constraint is always binding to check the accuracy of 
 * their linearization approach. Hopeully we can ind something that they missed. A main suspect 
 * is the asymmetry of policy functions.
 */

#define nk 512
#define nz 7
#define nxxi 7
#define nm1 1024 
#define tol 1e-6
#define maxiter 250
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

struct case1_hour {
	// Data Member are const coefficents and some model parameters
	double c0, c1, c_oneminusttheta, c_twominusttheta;
	para_struct para;

	// Construct a function of hour based on state and control variables 
	__host__ __device__
	case1_hour(double k, double z, double xxi, double m1, double zkttheta, para_struct _para) {
		c0 = (1-para.ddelta)*k*m1 - 1 + para.ddelta;
		c1 = (1-para.ddelta)*(1-k*m1-para.aalpha*para.ttheta/(1-para.ttheta));
		c_oneminusttheta = (1-1/xxi)*m1*zkttheta;
		c_twominusttheta = (1/xxi-1)*zkttheta*(m1+para.aalpha*para.ttheta/(k*(1-para.ttheta)));
		para = _para;
	};

	// The function of hour
	__host__ __device__
	double operator()(double n) {
		return c0 + c1*n + c_oneminusttheta*pow(n,1-para.ttheta) + c_twominusttheta*pow(n,2-para.ttheta);
	};

	__host__ __device__
	// The derivative of function
	double prime(double n) {
		return c1 + (1-para.ttheta)*c_oneminusttheta*pow(n,-para.ttheta) + (2-para.ttheta)*c_twominusttheta*pow(n,1-para.ttheta);
	};
};

struct case2_hour {
	// Data Member are const coefficents and some model parameters
	double c0,  c_oneminusttheta, c_minusttheta;
	para_struct para;

	// Construct a function of hour based on state and control variables 
	__host__ __device__
	case2_hour(double k, double z, double xxi, double m1, double zkttheta, para_struct _para) {
		double ddelta = para.ddelta;
		double ttheta = para.ttheta;
		double aalpha = para.aalpha;
		c0 = (1-ddelta)/m1;
		c_oneminusttheta = zkttheta*(ttheta/(m1*k)+(1-ttheta)/aalpha);
		c_minusttheta = -(1-ttheta)*zkttheta/aalpha;
		para = _para;
	};

	// The function of hour
	__host__ __device__
	double operator()(double n) {
		return c0 + c_oneminusttheta*pow(n,1-para.ttheta) + c_minusttheta*pow(n,-para.ttheta);
	};

	__host__ __device__
	// The derivative of function
	double prime(double n) {
		return (1-para.ttheta)*c_oneminusttheta*pow(n,-para.ttheta) + (-para.ttheta)*c_minusttheta*pow(n,-para.ttheta-1);
	};
};

// Define state struct that contain things known to shrink functor
struct state_struct {
	// Data member
	double k, z, xxi, m1, zkttheta;

	__host__ __device__
	void load(double _k, double _z, double _xxi, double _m1, double _zkttheta) {
		k = _k;
		z = _z;
		xxi = _xxi;
		m1 = _m1;
		zkttheta = _zkttheta;
	};
};



// Define a struct that contains variables implied by augmented state
struct control_struct {
	// Data member
	double kplus, c, n, w, d, mmu, Y, lhs1;

	// Constructor and finding the control variables
	__host__ __device__
	void compute(state_struct state, para_struct para, int binding) {
		if (binding == 1) {
			// Case 1: Binding
			double k = state.k;
			double z = state.z;
			double xxi = state.xxi;
			double m1 = state.m1;
			double zkttheta = state.zkttheta;
			n = newton(case1_hour(k,z,xxi,m1,zkttheta,para),0.0,1.0,0.3);
			Y = zkttheta*pow(n,1-para.ttheta);
			double MPK = para.ttheta*Y/k;
			kplus = Y/xxi;
			c = Y+(1-para.ddelta)*k-kplus;
			mmu = 1-(m1*c-1+para.ddelta)/MPK;	
			w = (1-mmu)*(1-para.ttheta)*Y/n;
			lhs1 = (1-mmu*xxi)/c;
			d = para.aalpha*c/(1-n);
		};

		if (binding == 0) {
			double k = state.k;
			double z = state.z;
			double xxi = state.xxi;
			double m1 = state.m1;
			double zkttheta = state.zkttheta;
			// Case 2: Not Binding
			n = newton(case2_hour(k,z,xxi,m1,zkttheta,para),0.0,1.0,0.3);
			Y = zkttheta*pow(n,1-para.ttheta);
			double MPK = para.ttheta*Y/k;
			mmu = 0;
			c = (1-para.ddelta+MPK)/m1;	
			kplus = (1-para.ddelta)*k + Y - c;
			w = (1-para.ttheta)*Y/n;
			lhs1 = (1-mmu*xxi)/c;
			d = para.aalpha*c/(1-n);
		};
	};
};

__host__ __device__
bool eureka(state_struct s, control_struct & u1, control_struct & u2, para_struct para, int i_z, int i_xxi, double* EM1_low, double* EM1_high, double* K) {
	double interp_low, interp_high;
	int i_kplus;

	// Case 1: Binding 
	u1.compute(s,para,1);

	// A series of tests whether it make senses
	if (u1.c <= 0) {
		goto case2;
	};
	if (u1.kplus < 0) {
		goto case2;
	};
	// if (u.kplus < K[0]) {
	// 	goto case2;
	// };
	// if (u.kplus > K[nk-1]) {
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
	interp_low = EM1_low[i_kplus+i_z*nk+i_xxi*nz*nk];
	interp_high = EM1_high[i_kplus+i_z*nk+i_xxi*nz*nk];
	if ( (u1.lhs1 > para.bbeta*interp_high) || (para.bbeta*interp_low > u1.lhs1) ) {
		// Euler equation 1 fails
		goto case2;
	};

	// Found a sensible m, break the loop
	return true;

case2: // Not Binding
	u2.compute(s,para,0);
	// A series of tests whether it make senses
	if (u2.c <= 0) {
		return false;
	};
	if (u2.kplus < 0) {
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
	interp_low = EM1_low[i_kplus+i_z*nk+i_xxi*nz*nk];
	interp_high = EM1_high[i_kplus+i_z*nk+i_xxi*nz*nk];
	if ( (u2.lhs1 > para.bbeta*interp_high) || (para.bbeta*interp_low > u2.lhs1) ) {
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
	para_struct para;

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
	para_struct _para)
	{
		K = K_ptr; Z = Z_ptr; XXI = XXI_ptr;
		V1_low = V1_low_ptr;
		V1_high = V1_high_ptr;
		Vplus1_low = Vplus1_low_ptr;
		Vplus1_high = Vplus1_high_ptr;
		EM1_low = EM1_low_ptr;
		EM1_high = EM1_high_ptr;
		flag = flag_ptr;
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

		// Find the "box" or "hypercube" that described m's range. Fancy word.
		double m1min = V1_low[index]; double m1max = V1_high[index];
		double m1min_old = m1min;
		double m1max_old = m1max;
		double step1 = (m1max-m1min)/double(nm1-1);
		double tempflag = 0.0;

		// Find and construct state and control, otherwise they won't update in the for loop
		double k =K[i_k]; double z=Z[i_z]; double xxi=XXI[i_xxi];
		double zkttheta = z*pow(k,para.ttheta); 
		state_struct s;
		control_struct u1, u2;

		// Initial search to find the min m
		for (int i_m1min = 0; i_m1min < nm1; i_m1min++) {
			// Construct state and find control variables
			double m1 = m1min+double(i_m1min)*step1;
			s.load(k,z,xxi,m1,zkttheta);
			if (eureka(s,u1,u2,para,i_z,i_xxi,EM1_low,EM1_high,K)) {
				m1min = m1;
				tempflag++;
				break;
			};
		};

		// Initial search to find the max m
		for (int i_m1max = 0; i_m1max < nm1; i_m1max++) {
			// Construct state and find control variables
			double m1 = m1max - double(i_m1max)*step1;
			s.load(k,z,xxi,m1,zkttheta);
			if (eureka(s,u1,u2,para,i_z,i_xxi,EM1_low,EM1_high,K)) {
				m1max = m1;
				tempflag++;
				break;
			};
		};

		// Update Vs
		flag[index] = double(tempflag)/double(nm1);
		if (tempflag == 0) {
			Vplus1_high[index] = m1max_old;
			Vplus1_low[index] = m1min_old; 
		} else {
			Vplus1_high[index] = m1max;
			Vplus1_low[index] = m1min;
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
	para.xxibar = 0.1;
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

	host_vector<double> h_V1_low(nk*nz*nxxi, 1/mkwidth*para.mkss);
	host_vector<double> h_V1_high(nk*nz*nxxi,mkwidth*para.mkss);
	host_vector<double> h_Vplus1_low(nk*nz*nxxi,1/mkwidth*para.mkss);
	host_vector<double> h_Vplus1_high(nk*nz*nxxi,mkwidth*para.mkss);

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
	double minK = 1/kwidth*para.kss;
	double maxK = kwidth*para.kss;
	linspace(9,11,nk,raw_pointer_cast(h_K.data()));

	// Create shocks grids
	host_vector<double> h_shockgrids(2*nz);
	double* h_shockgrids_ptr = raw_pointer_cast(h_shockgrids.data());
	double* h_P_ptr = raw_pointer_cast(h_P.data());
	gridgen_fptr linspace_fptr = &linspace; // select linspace as grid gen
	tauchen_vec(2,nz,5,para.A,para.Ssigma_e,h_shockgrids_ptr,h_P_ptr,linspace_fptr);
	for (int i_shock = 0; i_shock < nz; i_shock++) {
		h_Z[i_shock] = para.zbar*exp(h_shockgrids[i_shock+0*nz]);
		h_XXI[i_shock] = para.xxibar*exp(h_shockgrids[i_shock+1*nz]);
	};

	// Obtain initial guess from linear solution
	guess_linear(h_K, h_Z, h_XXI, h_V1_low, h_V1_high, para, 1/1.2, 1.2) ;

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
	while ((diff>tol)&&(iter<maxiter)&&(dist>0.005)){
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
				para)
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
	state_struct s;
	control_struct u;
	for (int i_k=0; i_k<nk; i_k++) {
		for (int i_z = 0; i_z < nz; i_z++) {
			for (int i_xxi=0; i_xxi < nxxi; i_xxi++) {
				int index = i_k+i_z*nk+i_xxi*nk*nz;
				double m1 = (h_V1_high[index]+h_V1_high[index])/2;
				double k = h_K[i_k];
				double z=h_Z[i_z]; double xxi=h_XXI[i_xxi];
				double zkttheta = z*pow(k,para.ttheta);
				s.load(k,z,xxi,m1,zkttheta);

				// Try not binding first
				u.compute(s,para,0);
				if (
						(s.xxi*u.kplus > u.Y) &&
						(u.c > 0) && 
						(u.kplus >= h_K[0]) &&
						(u.kplus <= h_K[nk-1]) &&
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
					u.compute(s,para,1);
					h_copt[index] = u.c;
					h_kopt[index] = u.kplus;
					h_nopt[index] = u.n;
					h_mmuopt[index] = u.mmu;
					h_dopt[index] = u.d;
					h_wopt[index] = u.w;
				};

				// Zoom in at the pike
				double step = (h_V1_high[index] - h_V1_low[index])/nm1;
				if ( (i_k==117) && (i_z==2) && (i_xxi==4)  ) {
					control_struct u1,u2;
					for (int i_m = 0; i_m < nm1; i_m++) {
						m1 = h_V1_low[index] + i_m*step;
						s.load(k,z,xxi,m1,zkttheta);
						u1.compute(s,para,1);
						u2.compute(s,para,0);
						h_kk_1[i_m] = u1.kplus;
						h_kk_2[i_m] = u2.kplus;
						h_lhs1_1[i_m] = u1.lhs1;
						h_lhs1_2[i_m] = u2.lhs1;

						// find euler related stuff
						int i_kplus_1 = fit2grid(u1.kplus,nk,h_K.data());
						int i_kplus_2 = fit2grid(u2.kplus,nk,h_K.data());
						h_rhslow_1[i_m] = para.bbeta*h_EM1_low[i_kplus_1+i_z*nk+i_xxi*nk*nz];
						h_rhshigh_1[i_m] = para.bbeta*h_EM1_high[i_kplus_1+i_z*nk+i_xxi*nk*nz];
						h_rhslow_2[i_m] = para.bbeta*h_EM1_low[i_kplus_2+i_z*nk+i_xxi*nk*nz];
						h_rhshigh_2[i_m] = para.bbeta*h_EM1_high[i_kplus_2+i_z*nk+i_xxi*nk*nz];
					};
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
	save_vec(h_lhs1_1,"./adrian_results/lhs1_1.csv");
	save_vec(h_lhs1_2,"./adrian_results/lhs1_2.csv");
	save_vec(h_rhslow_1,"./adrian_results/rhslow_1.csv");
	save_vec(h_rhshigh_1,"./adrian_results/rhshigh_1.csv");
	save_vec(h_rhslow_2,"./adrian_results/rhslow_2.csv");
	save_vec(h_rhshigh_2,"./adrian_results/rhshigh_2.csv");

	// Export parameters to MATLAB
	para.exportmatlab("./MATLAB/mypara.m");

	return 0;
}


