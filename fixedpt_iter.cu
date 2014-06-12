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

// Includes, Armadillo
#include "cppcode.h"

/* Includes, cuda */
#include <cublas_v2.h>
#include <cuda_helpers.h>

using namespace std;
using namespace thrust;

// #define M_PI 3.14159265358979323846264338328
#define nk 21
#define nb 1 
#define nz 11
#define nxxi 11 
#define nm1 501 
#define pk 7
#define pz 7
#define pxxi 7
#define tol 1e-10
#define maxiter 1e5
#define kwidth 0.3
#define bwidth 1.15 
#define mkwidth 15.0 // has to be a double 
#define llambda 0.5

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

void guess_linear(const host_vector<double> K, const host_vector<double> Z, const host_vector<double> XXI, host_vector<double> & M, para_struct para, double factor) {
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

	// Construct this object, create util from _util, etc.
	findnewM(double* K_ptr, double* K_cheby_ptr, double* Z_ptr, double* Z_cheby_ptr, double*XXI_ptr, double* XXI_cheby_ptr,
			 double* P_ptr, double* M_ptr, double* M_new_ptr, double* coeff_ptr, double _minK, double _maxK)
	{
		K = K_ptr; K_cheby = K_cheby_ptr;
		Z = Z_ptr; Z_cheby = Z_cheby_ptr;
		XXI = XXI_ptr; XXI_cheby = XXI_cheby_ptr;
		P = P_ptr;
		M = M_ptr; M_new = M_new_ptr;
		coeff = coeff_ptr;
		minK = _minK; maxK = _maxK;
	};

	__host__ __device__
	void operator()(int index) {
		int i_xxi = index/(nk*nz);
		int i_z = (index - i_xxi*nk*nz)/(nk);
		int i_k = index - i_xxi*nk*nz - i_z*nk;

		// Declare Variables
		double k, z, xxi, m;
		double kplus_cheby, zplus_cheby, xxiplus_cheby;	
		double d, mu, n, Y, c, kplus, w, ctilde;
		double EM;

		// Load Variables
		k = K[i_k]; z = Z[i_z]; xxi = XXI[i_xxi]; m = M[index];

		// Case 1: Not Binding
		mu = 0;
		n = newtonlabor(z,k,xxi,m,1);
		Y = z*pow(k,ttheta)*pow(n,1-ttheta);
		c = (1-ttheta)*Y/(aalpha*n);
		kplus = (1-ddelta)*k + Y - c;
		w = aalpha*c;
		d = c - w*n;
		kplus_cheby = -1 + (kplus-minK)/(maxK-minK)*(2);

		EM = 0;
		for (int i_zplus=0; i_zplus<nz; i_zplus++) {
			zplus_cheby = Z_cheby[i_zplus];
			for (int i_xxiplus=0; i_xxiplus<nxxi; i_xxiplus++) {
				xxiplus_cheby = XXI_cheby[i_xxiplus];
				EM += P[i_z+i_xxi*nz+nz*nxxi*i_zplus+nz*nxxi*nz*i_xxiplus]*chebyeval_multi(kplus_cheby,zplus_cheby,xxiplus_cheby,coeff);
			};
		};
		ctilde = (1-mu*xxi)/(bbeta*EM);

		// Check whether implied policy functions make sense
		if (
			(c>0) && (ctilde>0) && (xxi*kplus>w*n) && (kplus > minK) && (kplus<maxK) && (n>0) && (n<=1)
		   )
		{
			M_new[index] = (z*ttheta*pow(k,ttheta-1)*pow(n,1-ttheta)+1-ddelta)/ctilde;
			goto stop;
		} else {
			goto case2;
		}; 

		// Case 2: Binding
		case2:
		n = newtonlabor(z,k,xxi,m,2); 
		c = (1-ddelta)/m + ttheta*z*pow(k,ttheta)*pow(n,1-ttheta)/(m*k); // correct
		w = aalpha*c;	// correct
		kplus = w*n/xxi;	// correct
		d = c - w*n;	// correct
		mu = (1-ttheta)*z*pow(k,1-ttheta)*pow(n,-ttheta)/w - 1;	// correct
		kplus_cheby = -1 + (kplus-minK)/(maxK-minK)*(2);

		EM = 0;
		for (int i_zplus=0; i_zplus<nz; i_zplus++) {
			zplus_cheby = Z_cheby[i_zplus];
			for (int i_xxiplus=0; i_xxiplus<nxxi; i_xxiplus++) {
				xxiplus_cheby = XXI_cheby[i_xxiplus];
				EM += P[i_z+i_xxi*nz+nz*nxxi*i_zplus+nz*nxxi*nz*i_xxiplus]*chebyeval_multi(kplus_cheby,zplus_cheby,xxiplus_cheby,coeff);
			};
		};
		ctilde = (1-mu*xxi)/(bbeta*EM);

		if (
			(ctilde>0) && (c>0) && (kplus>minK) && (kplus<maxK) && (mu>=0) && (n>0) && (n<=1)
		   )
		{
			M_new[index] = (z*ttheta*pow(k,ttheta-1)*pow(n,1-ttheta)+1-ddelta)/ctilde;
		} else {
			printf("No solution at k=%f, z=%f, xxi=%f, m=%f.\nPolicies are: c=%f, ctilde=%f, kplus=%f, mu=%f, n=%f\n====================================================\n",k,z,xxi,m,c,ctilde,kplus,mu,n);
			// M_new[index] = NAN;
		};

		// The end here
		stop:
		int shit = 0;
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
int main()
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
	host_vector<double> h_K_cheby(nk); 
	host_vector<double> h_Z(nz);
	host_vector<double> h_Z_cheby(nz);
	host_vector<double> h_XXI(nxxi);
	host_vector<double> h_XXI_cheby(nxxi);
	host_vector<double> h_P(nz*nxxi*nz*nxxi);
	host_vector<double> h_flag(nk*nb*nz*nxxi, 0); 
	host_vector<double> h_X(nk*nz*nxxi*(1+pk)*(1+pz)*(1+pxxi)); 
	host_vector<double> h_projector((1+pk)*(1+pz)*(1+pxxi)*nk*nz*nxxi); 
	host_vector<double> h_M(nk*nz*nxxi); 
	host_vector<double> h_M_new(nk*nz*nxxi); 
	host_vector<double> h_coeff((1+pk)*(1+pz)*(1+pxxi),0.1); 
	host_vector<double> h_coeff_temp((1+pk)*(1+pz)*(1+pxxi),0.1); 
	host_vector<double> h_coeff_new((1+pk)*(1+pz)*(1+pxxi),0.1); 
	
	// Create capital grid
	double* h_K_cheby_ptr = raw_pointer_cast(h_K_cheby.data());
	chebyroots(nk,h_K_cheby_ptr);
	h_K = h_K_cheby;
	double* h_K_ptr = raw_pointer_cast(h_K.data());
	double minK = (1-kwidth)*para.kss;
	double maxK = (1+kwidth)*para.kss;
	cout << "minK: " << minK << endl;
	cout << "maxK: " << maxK << endl;
	fromchebydomain(minK, maxK, nk, h_K_ptr);

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
	double* h_Z_cheby_ptr = raw_pointer_cast(h_Z_cheby.data());
	double* h_XXI_cheby_ptr = raw_pointer_cast(h_XXI_cheby.data());
	chebyroots(nz,h_Z_cheby_ptr);
	chebyroots(nxxi,h_XXI_cheby_ptr);

	// Create Initial M generated from linear solution
	guess_linear(h_K, h_Z, h_XXI, h_M, para, 1.0) ;

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
	double diff = 10; double dist; int iter = 0;
	while ((diff>tol)&&(iter<maxiter)){
		// Find the current M at each state. Does y = A*x
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
				findnewM(d_K_ptr, d_K_cheby_ptr, d_Z_ptr, d_Z_cheby_ptr, d_XXI_ptr, d_XXI_cheby_ptr,
					     d_P_ptr, d_M_ptr, d_M_new_ptr, d_coeff_ptr, minK, maxK)
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
		printf("=======================================================\n=== Iteration No. %i finished \n=======================================================\n",
				iter);
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
	h_coeff = d_coeff;
	ofstream fout_coeff("coeff.csv", ios::trunc);
	for (int index=0; index<(1+pk)*(1+pz)*(1+pxxi); index++) {
		fout_coeff << setprecision(18) << h_coeff[index] << '\n';
	};
	// ofstream fout_V1_low("V1_low.csv", ios::trunc); ofstream fout_V1_high("V1_high.csv", ios::trunc);
	// ofstream fout_kopt("koptcuda.csv", ios::trunc); ofstream fout_copt("coptcuda.csv", ios::trunc);
	// ofstream fout_R("Rcuda.csv", ios::trunc);
	// ofstream fout_d("dcuda.csv", ios::trunc); ofstream fout_n("ncuda.csv", ios::trunc);
	// ofstream fout_wage("wagecuda.csv", ios::trunc);
	// ofstream fout_Kgrid("Kgridcuda.csv", ios::trunc);
	// ofstream fout_Zgrid("Zgridcuda.csv", ios::trunc); ofstream fout_XIgrid("XIgridcuda.csv", ios::trunc);
	// ofstream fout_mmu("mmucuda.csv", ios::trunc); 
	// ofstream fout_P("Pcuda.csv", ios::trunc);
	// ofstream fout_flag("flagcuda.csv", ios::trunc);
	// 
	// for (int index=0; index<nk*nz*nxi; index++) {
	// 	fout_V1_low << setprecision(16) << h_V1_low[index] << '\n';
	// 	fout_V1_high << setprecision(16) << h_V1_high[index] << '\n';
	// 	fout_flag << setprecision(16) << h_flag[index] << '\n';
	// 	int i_xi = index/(nk*nz);
	// 	int i_z  = (index-i_xi*nk*nz)/(nk);
	// 	int i_k = index - i_xi*nk*nz - i_z*nk ;
	// 	double m1_high = (h_V1_high[index]+h_V1_high[index])/2;
	// 	double m1_low = (h_V1_low[index]+h_V1_low[index])/2;
	// 	double m1 = m1_high;
	// 	double k =h_K[i_k]; 
	// 	double z=h_Z[i_z]; double xi=h_XI[i_xi];
	// 	

	// 	double mu; double w; double c; double kplus; double d; double n;
	// 	double Y;


	// 	// Case 2: Binding
	// 	n = newtonlabor(z,k,xi,m1,2); 
	// 	c = (1-ddelta)/m1 + z*ttheta*pow(k,ttheta-1)*pow(n,1-ttheta)/m1; // correct
	// 	w = aalpha*c;	// correct
	// 	kplus = w*n/xi;	// correct
	// 	d = c - w*n;	// correct
	// 	mu = (1-ttheta)*z*pow(k/n,ttheta)/w - 1;	// correct
	// 	if (mu<0) {
	// 		mu = 0;
	// 		n = newtonlabor(z,k,xi,m1,1);
	// 		Y = z*pow(k,ttheta)*pow(n,1-ttheta);
	// 		c = (1-ttheta)*Y/(aalpha*n);
	// 		kplus = (1-ddelta)*k + Y - c;
	// 		w = aalpha*c;
	// 		d = c - w*n;
	// 	};

	// 	fout_kopt << setprecision(16) << kplus << '\n';
	// 	fout_copt << setprecision(16) <<  c << '\n';
	// 	fout_R << setprecision(16) <<  1 << '\n';
	// 	fout_d << setprecision(16) <<  d << '\n';
	// 	fout_n << setprecision(16) <<  n << '\n';
	// 	fout_wage << setprecision(16) <<  w << '\n';
	// 	fout_mmu << setprecision(16) <<  mu << '\n';
	// };
	// 
	// for (int i=0; i<nk; i++) {
	// 	fout_Kgrid << setprecision(16) << h_K[i] << '\n';
	// };
	// for (int i=0; i<nz; i++) {
	// 	fout_Zgrid << setprecision(16) << h_Z[i] << '\n';
	// };
	// for (int i=0; i<nxi; i++) {
	// 	fout_XIgrid << setprecision(16) << h_XI[i] << '\n';
	// };
	// for (int i=0; i<nz*nxi*nz*nxi; i++) {
	// 	fout_P << setprecision(16) << h_P[i] << '\n';
	// };

	// fout_V1_low.close(); fout_V1_high.close();
	// fout_Kgrid.close();
	// fout_Zgrid.close(); fout_XIgrid.close();

	// //  fout_kopt_high.close(); fout_copt_high.close();
	// //  fout_R_high.close(); fout_d_high.close();
	// //  fout_n_high.close(); fout_mmu_high.close();

	// fout_kopt.close(); fout_copt.close();
	// fout_R.close(); fout_d.close();
	// fout_n.close(); fout_mmu.close();
	// fout_wage.close();

	// fout_P.close();
	// fout_flag.close();

	// printstuff< host_vector<double> >(d_X);
	return 0;
};


