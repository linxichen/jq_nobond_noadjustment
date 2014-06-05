/* This is a CUDA implementation of Jermann Quadrini 2013 AER  
 * We simulate a equilibrium that collateral constraint is always binding to check the accuracy of 
 * their linearization approach. Hopeully we can ind something that they missed. A main suspect 
 * is the asymmetry of policy functions.
 */

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

using namespace std;
using namespace thrust;

// Define an class that contains parameters and steady states
struct para_struct {
	// Accuracy controls
	int nk ;
	int nb ;
	int nz ;
	int nxxi;
	int nm1;
	int maxiter;
	double tol;
	double kwidth ;
	double bwidth ;
	double mkwidth;

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
		fileout << setprecision(16) << "nb=" << nb << ";"<< endl;
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
	for (int i_k=0; i_k<para.nk; i_k++) {
		for (int i_z = 0; i_z < para.nz; i_z++) {
			for (int i_xxi = 0; i_xxi < para.nxxi; i_xxi++) {
				double temp = para.mkss+Pphi[8+0*9]*(K[i_k]-para.kss) + Pphi[8+1*9]*(log(Z[i_z])-log(para.zbar))+ Pphi[8+2*9]*(log(XXI[i_xxi])-log(para.xxibar));
				V1_low[i_k+para.nk*i_z+para.nk*para.nz*i_xxi] = factor_low*temp;
				V1_high[i_k+para.nk*i_z+para.nk*para.nz*i_xxi] = factor_high*temp;
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
		c1 = (1-para.ddelta)*(1-k*m1+para.aalpha*para.ttheta/(1-para.ttheta));
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
		c_oneminusttheta = (ttheta*zkttheta/(m1*k)+(1-ttheta)*zkttheta/aalpha);
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
		return c0 + (1-para.ttheta)*c_oneminusttheta*pow(n,-para.ttheta) + (-para.ttheta)*c_minusttheta*pow(n,-para.ttheta-1);
	};
};

// Eureka function check whether a tuple (STATE,SHOCK,SHADOW) can survive to next iteration
__host__ __device__
bool eureka(double k, double z, double xxi,
            double m1, double zkttheta, int i_z, int i_xi,
            double* K, 
			double* EM1_low, double* EM1_high, 
			para_struct para) {
	
	// Declare control variables
	double n, Y, MPK, kplus, c, mmu, w, lhs1;
	int i_kplus;
	double interp_low, interp_high;

	// Case 1: Binding
	n = newton(case1_hour(k,z,xxi,m1,zkttheta,para),0.0,1.0,0.3);
	// printf("Hour solved here is: %f\n",n);
	Y = zkttheta*pow(n,1-para.ttheta);
	MPK = para.ttheta*Y/k;
	kplus = Y/xxi;
	c = Y+(1-para.ddelta)*k-kplus;
	mmu = 1-(m1*c-1+para.ddelta)/MPK;	
	w = (1-mmu)*(1-para.ttheta)*Y/n;
	// d = c-w*n;
	i_kplus = fit2grid(kplus,para.nk,K);
	lhs1 = (1-xxi*mmu)/c;
	// interp_low = EM1_low[i_kplus+para.nk*(i_z+i_xi*para.nz)] + (kplus-K[i_kplus])*(EM1_low[i_kplus+1+para.nk*(i_z+i_xi*para.nz)]-EM1_low[i_kplus+para.nk*(i_z+i_xi*para.nz)])/(K[i_kplus+1]-K[i_kplus]);
	// interp_high = EM1_high[i_kplus+para.nk*(i_z+i_xi*para.nz)] + (kplus-K[i_kplus])*(EM1_high[i_kplus+1+para.nk*(i_z+i_xi*para.nz)]-EM1_high[i_kplus+para.nk*(i_z+i_xi*para.nz)])/(K[i_kplus+1]-K[i_kplus]);
	interp_low = EM1_low[i_kplus+i_z*para.nk+i_xi*para.nz*para.nk];
	interp_high = EM1_high[i_kplus+i_z*para.nk+i_xi*para.nz*para.nk];
	
	if (
		(para.bbeta*interp_low <= lhs1) &&
		(lhs1 <=para.bbeta*interp_high) &&
		(c>0) && (mmu>=0) && (w>=0) && (n>0) && (n<1)//  &&  
	   )
	{
		return true;
	};

	// Case 2: Not Binding
	n = newton(case2_hour(k,z,xxi,m1,zkttheta,para),0.0,1.0,0.3);
	Y = zkttheta*pow(n,1-para.ttheta);
	MPK = para.ttheta*Y/k;
	mmu = 0;
	c = (1-para.ddelta+MPK)/m1;	
	kplus = (1-para.ddelta)*k + Y - c;
	w = (1-para.ttheta)*Y/n;
	lhs1 = (1-xxi*mmu)/c;
	i_kplus = fit2grid(kplus,para.nk,K);
	// interp_low = EM1_low[i_kplus+para.nk*(i_z+i_xi*para.nz)] + (kplus-K[i_kplus])*(EM1_low[i_kplus+1+para.nk*(i_z+i_xi*para.nz)]-EM1_low[i_kplus+para.nk*(i_z+i_xi*para.nz)])/(K[i_kplus+1]-K[i_kplus]);
	// interp_high = EM1_high[i_kplus+para.nk*(i_z+i_xi*para.nz)] + (kplus-K[i_kplus])*(EM1_high[i_kplus+1+para.nk*(i_z+i_xi*para.nz)]-EM1_high[i_kplus+para.nk*(i_z+i_xi*para.nz)])/(K[i_kplus+1]-K[i_kplus]);
	interp_low = EM1_low[i_kplus+i_z*para.nk+i_xi*para.nz*para.nk];
	interp_high = EM1_high[i_kplus+i_z*para.nk+i_xi*para.nz*para.nk];
	if (
		(para.bbeta*interp_low <= lhs1) &&
		(lhs1 <=para.bbeta*interp_high) &&
		(c>0) && (w>=0) && (xxi*kplus>Y) && (n>0) && (n<1)
	   )
	{
		return true;
	} else {
		return false;
	};
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
		int subs[3];
		int size_vec[3];
		size_vec[0] = para.nk;
		size_vec[1] = para.nz;
		size_vec[2] = para.nxxi;
		ind2sub(3,size_vec,index,subs);
		int i_k = subs[0];
		int i_z = subs[1];
		int i_xxi = subs[2];
		double k =K[i_k]; double z=Z[i_z]; double xxi=XXI[i_xxi];
		double zkttheta = z*pow(k,para.ttheta); 

		// Find the "box" or "hypercube" that described m's range. Fancy word.
		double m1min = V1_low[index]; double m1max = V1_high[index];
		double m1min_old = m1min;
		double m1max_old = m1max;
		double step1 = (m1max-m1min)/double(para.nm1-1);
		double tempflag = 0.0;
		int nm1 = para.nm1;

		// Initial search to find the min m
		for (int m_index = 0; m_index < nm1; m_index++) {
			int i_m1 = m_index;
			double m1=m1min+double(i_m1)*step1;
			if (eureka(k,z,xxi,m1,zkttheta,i_z,i_xxi,K,
				EM1_low,EM1_high,para))
			{
				tempflag++;
				m1min = m1-step1;
				break; // break and only break one layer of for loop
			};
		};

		// // "Trace-back" to refine the min_m1, assuming we found at least one m !!!
		// double min_step = step1/(nm1-1);
		// double m1min_left = m1min - step1; 
		// for (int i_m1min = 0; i_m1min < nm1; i_m1min++) {
		// 	double m1 = m1min_left + i_m1min*min_step;
		// 	if (eureka(k,z,xxi,m1,zkttheta,i_z,i_xxi,K,
		// 		EM1_low,EM1_high,para))
		// 	{
		// 		tempflag++;
		// 		m1min = m1;
		// 		break; // break and only break one layer of for loop
		// 	};
		// };

		// Initial search to find the max m
		for (int m_index = nm1-1; m_index >= 0; m_index--) {
			int i_m1 = m_index;
			double m1=m1min+double(i_m1)*step1;
			if (eureka(k,z,xxi,m1,zkttheta,i_z,i_xxi,K,
				EM1_low,EM1_high,para))
			{
				tempflag++;
				m1max = m1+step1;
				break; // break and only break one layer of for loop
			};
		};

		// // "Trace-back" to refine the max_m1, assuming we found at least one m !!!
		// double max_step = step1/(nm1-1);
		// double m1max_right = m1max + step1; 
		// for (int i_m1max = 0; i_m1max < nm1; i_m1max++) {
		// 	double m1 = m1max_right - i_m1max*max_step;
		// 	if (eureka(k,z,xxi,m1,zkttheta,i_z,i_xxi,K,
		// 		EM1_low,EM1_high,para))
		// 	{
		// 		tempflag++;
		// 		m1max = m1;
		// 		break; // break and only break one layer of for loop
		// 	};
		// };

		// Update Vs
		flag[index] = double(tempflag)/double(nm1);
		if (tempflag == 0) {
			Vplus1_high[index] = m1min_old;
			Vplus1_low[index] = m1max_old; 
		} else {
			Vplus1_high[index] = m1max;
			Vplus1_low[index] = m1min;
		};
	}
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

int main()
{
	// Initialize Parameters
	para_struct para;

	// Set Accuracy Parameters
	para.nk = 256;
	para.nb = 1 ;
	para.nz = 7;
	para.nxxi = 7;
	para.nm1 = 2560 ;
	para.tol = 1e-6;
	para.maxiter = 1e5;
	para.kwidth = 1.5 ;
	para.bwidth = 1.15 ;
	para.mkwidth = 20.0 ; 

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
	cout << setprecision(16) << "tol: " << para.tol << endl;

	// Select Device
	int num_devices;
	cudaGetDeviceCount(&num_devices);
	if (num_devices>1) {
		cudaSetDevice(1);
	};
	// Only for cuBLAS
	const double alpha = 1.0;
	const double beta = 0.0;

	// Create all STATE, SHOCK grids here
	host_vector<double> h_K(para.nk); 
	host_vector<double> h_Z(para.nz);
	host_vector<double> h_XXI(para.nxxi);

	host_vector<double> h_V1_low(para.nk*para.nz*para.nxxi, 1/para.mkwidth*para.mkss);
	host_vector<double> h_V1_high(para.nk*para.nz*para.nxxi,para.mkwidth*para.mkss);
	host_vector<double> h_Vplus1_low(para.nk*para.nz*para.nxxi,1/para.mkwidth*para.mkss);
	host_vector<double> h_Vplus1_high(para.nk*para.nz*para.nxxi,para.mkwidth*para.mkss);

	host_vector<double> h_EM1_low(para.nk*para.nz*para.nxxi,0.0);
	host_vector<double> h_EM1_high(para.nk*para.nz*para.nxxi,0.0);

	host_vector<double> h_P(para.nz*para.nxxi*para.nz*para.nxxi, 1.0/double(para.nz*para.nxxi));
	host_vector<double> h_flag(para.nk*para.nz*para.nxxi, 0); 

	// host_vector<double> h_c(para.nk*para.nz*para.nxxi*para.nm1);
	// host_vector<double> h_n(para.nk*para.nz*para.nxxi*para.nm1);
	// host_vector<double> h_kplus(para.nk*para.nz*para.nxxi*para.nm1);
	// host_vector<double> h_mmu(para.nk*para.nz*para.nxxi*para.nm1);
	// host_vector<double> h_lhs(para.nk*para.nz*para.nxxi*para.nm1);
	// host_vector<double> h_rhs_low(para.nk*para.nz*para.nxxi*para.nm1);
	// host_vector<double> h_rhs_high(para.nk*para.nz*para.nxxi*para.nm1);
	
	// Create capital grid
	double minK = 1/para.kwidth*para.kss;
	double maxK = para.kwidth*para.kss;
	linspace(minK,maxK,para.nk,raw_pointer_cast(h_K.data()));
	save_vec(h_K,para.nk,"Kgrid.csv");

	// Create shocks grids
	host_vector<double> h_shockgrids(2*para.nz);
	double* h_shockgrids_ptr = raw_pointer_cast(h_shockgrids.data());
	double* h_P_ptr = raw_pointer_cast(h_P.data());
	gridgen_fptr linspace_fptr = &linspace; // select linspace as grid gen
	tauchen_vec(2,para.nz,3,para.A,para.Ssigma_e,h_shockgrids_ptr,h_P_ptr,linspace_fptr);
	for (int i_shock = 0; i_shock < para.nz; i_shock++) {
		h_Z[i_shock] = para.zbar*exp(h_shockgrids[i_shock+0*para.nz]);
		h_XXI[i_shock] = para.xxibar*exp(h_shockgrids[i_shock+1*para.nz]);
	};
	display_vec(h_Z);
	display_vec(h_XXI);
	save_vec(h_Z,"Zgrid.csv");
	save_vec(h_XXI,"XXIgrid.csv");
	save_vec(h_P,"Pcuda.csv");

	// Obtain initial guess from linear solution
	guess_linear(h_K, h_Z, h_XXI, h_V1_low, h_V1_high, para, 0.5, 1.5) ;
	save_vec(h_V1_low,"V1_low_guess.csv");
	save_vec(h_V1_high,"V1_high_guess.csv");

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
	counting_iterator<int> end(para.nk*para.nz*para.nxxi);

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
	
	double diff = 10; double dist; int iter = 0;
	while ((diff>para.tol)&&(iter<para.maxiter)){
		// Find EMs for low and high 
		cublasDgemm(handle,
			CUBLAS_OP_N,  
			CUBLAS_OP_T,
			para.nk*para.nb, para.nz*para.nxxi, para.nz*para.nxxi,
			&alpha,
			d_V1_low_ptr, 
			para.nk*para.nb, 
			d_P_ptr,
			para.nz*para.nxxi,
			&beta,
			d_EM1_low_ptr,
			para.nk*para.nb);
		cublasDgemm(handle,
			CUBLAS_OP_N,  
			CUBLAS_OP_T,
			para.nk*para.nb, para.nz*para.nxxi, para.nz*para.nxxi,
			&alpha,
			d_V1_high_ptr, 
			para.nk*para.nb, 
			d_P_ptr,
			para.nz*para.nxxi,
			&beta,
			d_EM1_high_ptr,
			para.nk*para.nb);

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
		cout << "Vplus1[100-1] (the spike) range is " << d_Vplus1_low[100-1] << ", " << d_Vplus1_high[100-1] << endl;

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
	h_flag = d_flag;
	save_vec(h_V1_low,"V1_low.csv");
	save_vec(h_V1_high,"V1_high.csv");
	save_vec(h_flag,"flagcuda.csv");
	
	ofstream fout_kopt("koptcuda.csv", ios::trunc); ofstream fout_copt("coptcuda.csv", ios::trunc);
	ofstream fout_R("Rcuda.csv", ios::trunc);
	ofstream fout_wage("wagecuda.csv", ios::trunc);
	ofstream fout_d("dcuda.csv", ios::trunc); ofstream fout_n("ncuda.csv", ios::trunc);
	ofstream fout_Kgrid("Kgridcuda.csv", ios::trunc);
	ofstream fout_Zgrid("Zgridcuda.csv", ios::trunc); ofstream fout_XXIgrid("XXIgridcuda.csv", ios::trunc);
	ofstream fout_mmu("mmucuda.csv", ios::trunc); ofstream fout_P("Pcuda.csv", ios::trunc);
	
	double ttheta = para.ttheta;
	int nk = para.nk;
	int nz = para.nz;
	int nxxi = para.nxxi;
	for (int index=0; index<nk*nz*nxxi; index++) {
		int i_xxi = index/(nk*nz);
		int i_z  = (index-i_xxi*nk*nz)/(nk);
		int i_k = index - i_xxi*nk*nz - i_z*nk ;
		double m1 = (h_V1_low[index]+h_V1_low[index])/2;
		double k =h_K[i_k];
		double z=h_Z[i_z]; double xxi=h_XXI[i_xxi];
		double zkttheta = z*pow(k,ttheta);

		// Declare control variables
		double n, Y, MPK, kplus, c, mmu, w, d;

		// Case 1: Binding
		n = newton(case1_hour(k,z,xxi,m1,zkttheta,para),0.0,1.0,0.3);
		Y = zkttheta*pow(n,1-para.ttheta);
		MPK = para.ttheta*Y/k;
		kplus = Y/xxi;
		c = Y+(1-para.ddelta)*k-kplus;
		mmu = 1-(m1*c-1+para.ddelta)/MPK;	
		w = (1-mmu)*(1-para.ttheta)*Y/n;
		d = c - w*n;
	
		if (mmu<0)
		{
			// Case 2: Not Binding
			n = newton(case2_hour(k,z,xxi,m1,zkttheta,para),0.0,1.0,0.3);
			Y = zkttheta*pow(n,1-para.ttheta);
			MPK = para.ttheta*Y/k;
			mmu = 0;
			c = (1-para.ddelta+MPK)/m1;	
			kplus = (1-para.ddelta)*k + Y - c;
			w = (1-para.ttheta)*Y/n;
			d = c - w*n;
		};

		fout_kopt << kplus << '\n';
		fout_copt << c << '\n';
		fout_R << 1 << '\n';
		fout_d << d << '\n';
		fout_n << n << '\n';
		fout_mmu << mmu << '\n';
		fout_wage << w << '\n';
	};
	
	for (int i=0; i<nk; i++) {
		fout_Kgrid << h_K[i] << '\n';
	};
	for (int i=0; i<nz; i++) {
		fout_Zgrid << h_Z[i] << '\n';
	};
	for (int i=0; i<nxxi; i++) {
		fout_XXIgrid << h_XXI[i] << '\n';
	};
	for (int i=0; i<nz*nxxi*nz*nxxi; i++) {
		fout_P << h_P[i] << '\n';
	};

	fout_Kgrid.close();
	fout_Zgrid.close(); fout_XXIgrid.close();
	fout_kopt.close(); fout_copt.close();
	fout_R.close(); fout_d.close();
	fout_n.close(); fout_mmu.close(); fout_P.close();
	fout_wage.close();

	// Export parameters to MATLAB
	para.exportmatlab("./MATLAB/mypara.m");
	return 0;
}


