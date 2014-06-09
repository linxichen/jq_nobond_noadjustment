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

// Parameters from the AER paper
#define aalpha 2.7189
#define bbeta 0.9825
#define ddelta 0.025
#define ttheta 0.36
// #define kkappa 0.1460
#define kkappa 0.0265
// #define ttau 0.3500
// #define xxibar 0.1634
// #define xxibar 0.1634
// #define xxibar 0.02695
#define xxibar 0.085
#define zbar 1.0
// #define dbar 2.56514
#define uupsilon 1.0
#define oomega 1.4
// #define aalpha 1.3659
#define rrhozz 0.8147
#define rrhozxxi -0.1234
#define rrhoxxiz 0.15
#define rrhoxxixxi 0.8393
#define ssigmaz 0.0042
#define ssigmaxxi 0.0072

// Evaluate Chebychev polynomial of any degree
__host__ __device__
double chebypoly(const int p, const double x) {
	switch (p) {
		case 0: // 0-th order Chebyshev Polynomial
			return 1;
		case 1:
			return x;
		case 2:
			return 2*x*x - 1;
		case 3:
			return 4*x*x*x - 3*x;
	}
	
	// When p>=4, apply the recurrence relation
	double lag1 = 4*x*x*x -3*x;
	double lag2 = 2*x*x - 1;
	double lag0;
	int distance = p - 3;
	while (distance >= 1) {
		lag0 = 2*x*lag1 - lag2;
		lag2 = lag1;
		lag1 = lag0;
		distance--;
	};
	return lag0;
};

// Evaluate Chebychev polynomial of any degree
__host__ __device__
int chebyroots_cuda(const int p, double* roots) {
	for (int i=0; i<p; i++) {
		double stuff = p - 0.5 - 1*i;
		roots[i] = cos(M_PI*(stuff)/(p));
	};

	// Account for the fact that cos(pi/2) is not exactly zeros
	if (p%2) {
		roots[(p-1)/2] = 0;
	};
	return 0;
};

void ind2sub_cuda(int length_size, device_vector<int> siz_vec, int index, device_vector<int> subs) {
// Purpose:     Converts index to subscripts. i -> [i_1, i_2, ..., i_n]
//
// Input:       length_size = # of coordinates, i.e. how many subscripts you are getting
//              siz_vec = vector that stores the largest coordinate value for each subscripts. Or the dimensions of matrices
//              index = the scalar index
//
// Ouput:       subs = the vector stores subscripts
    device_vector<int> cumdim(length_size);
    cumdim[0] = 1;
    for (int i=1; i<length_size; i++) {
        cumdim[i] = cumdim[i-1]*siz_vec[i-1];
    };
    int done = 0;
    for ( int i=length_size-1; i>=0; i-- ) {
        subs[i] = (index - done)/cumdim[i];
        done += subs[i]*cumdim[i];
    };
};


// Evaluate Chebychev approximation of any degree
__host__ __device__
double chebyeval(int p, double x, double* coeff) {
	// Note that coefficient vector has p+1 values
	double sum = 0;
	for (int i=0; i<=p; i++) {
		sum += coeff[i]*chebypoly(i,x);	
	};
	return sum;
};

// This is a very model-specific version of the function it should be.
// Future modification needs to deal with the creation of temporary array somehow.
__host__ __device__
double chebyeval_multi ( double k_cheby, double z_cheby, double xxi_cheby, double* coeff ) {
	double eval = 0;
	for (int t_xxi=0; t_xxi <= pxxi; t_xxi++) {
		for (int t_z=0; t_z <= pz; t_z++) {
			for (int t_k=0; t_k <=pk; t_k++ ) {
				eval += coeff[t_k+t_z*(1+pk)+t_xxi*(1+pk)*(1+pz)]*
				        chebypoly(t_k,k_cheby)*chebypoly(t_z,z_cheby)*chebypoly(t_xxi,xxi_cheby);
			};
		};
	};
	return eval;
};

// __host__ __device__
// double chebyeval_multi(device_vector<int> p_vec, device_vector<double> x_vec, device_vector<double> coeff) {
// 	int n_var = p_vec.size();
// 	device_vector<int> size_vec(n_var);
// 	device_vector<int> subs(n_var);
// 	int tot_terms = 1;
// 	// Do this because for each variables there's p+1 coefficients!
// 	for (int i_var = 0; i_var < n_var; i_var++) {
// 		tot_terms *= (p_vec[i_var]+1);
// 	};
// 
// 	// 
// 	double sum = 0;
// 	for (int index = 0; index < tot_terms; index++) {
// 		ind2sub_cuda(n_var, size_vec, index, subs);
// 		printf("%d",subs[0]);
// 	};
// 
// };

// This function fit a valuex x to a grid X of size n
__host__ __device__
int fit2grid(const double x, const int n, const double* X) {
	if (x < X[0]) {
		return 0;
	} else if (x > X[n-1]) {
		return n-1;
	} else {
		int left=0; int right=n-1; int mid=(n-1)/2;
		while(right-left>1) {
			mid = (left + right)/2;
			if (X[mid]==x) {
				return mid;
			} else if (X[mid]<x) {
				left = mid;
			} else {
				right = mid;
			};

		};
		return left;
	}
};

// The objective function for case 2
__host__ __device__
double objective2(double z,double k, double xxi, double m1, double n) {
	return (1-ddelta)*k-(1-ddelta)/m1-(1-ddelta)*aalpha*n/(xxi*m1)
		+ z*pow(k,ttheta)*pow(n,1-ttheta)*(1-ttheta/(m1*k)-ttheta*aalpha*n/(m1*xxi*k));
};

// The  derivative for case 2
__host__ __device__
double derivative2(double z,double k,double xxi,double m1,double n) {
	return -(1-ddelta)*aalpha/(xxi*m1)+(1-ttheta)*z*pow(k,ttheta)*pow(n,-ttheta)*
		(1-ttheta/(k*m1)-ttheta*aalpha*n/(m1*xxi*k)) + z*pow(k,ttheta)*pow(n,1-ttheta)*(-ttheta*aalpha/(m1*xxi*k));
};

// The objective function for case 1
__host__ __device__
double objective1(double z,double k,double xxi,double m1,double n) {
	return z*pow(k,ttheta)*pow(n,-ttheta)*(m1*(1-ttheta)/aalpha -ttheta*n/k) -1 + ddelta;
};

// The  derivative for case 1
__host__ __device__
double derivative1(double z,double k,double xxi,double m1,double n) {
	return z*pow(k,ttheta)*pow(n,-ttheta-1)*(-ttheta*(1-ttheta)*m1/aalpha+(ttheta*ttheta-ttheta)*n/k);
};

// Apply Newton's Method to find labor hours 
__host__ __device__
double newtonlabor(double z,double k,double xxi,double m1,int lever) {
	double n_old = 0.3;
	double n_new = 0.3;
	double mytol = 0.0000001;
	int mymaxiter = 30;
	int found = 0; // 1 means we found the solution

	if (lever==1) {
		// Check if it's possible to find a root
		double f_old = objective1(z,k,xxi,m1,0); 
		double f_new = objective1(z,k,xxi,m1,1); 

		if (f_old*f_new>0) {
			return -99.99;
		};

		int iter = 0;
		while ( (found != 1) && (iter <= mymaxiter) ) {
			n_new = n_old - objective1(z,k,xxi,m1,n_old)/derivative1(z,k,xxi,m1,n_old);
			if (n_new<0) {
				n_new = 0.0;
			};
			if (n_new>1) {
				n_new = 1.0;
			};
			double error = abs((n_new-n_old)/n_new);
			if (error<mytol) {
				found = 1;
			} else {
				n_old = n_new;
			};

			iter++;
		};

		if (found==1) {
			return n_new;
		} else {
			return -99.99;
		};
	};

	// Check if it's possible to find a root
	double f_old = objective2(z,k,xxi,m1,0); 
	double f_new = objective2(z,k,xxi,m1,1); 

	if (f_old*f_new>0) {
		return -99.99;
	};

	int iter = 0;
	while ( (found != 1) && (iter <= mymaxiter) ) {
		n_new = n_old - objective2(z,k,xxi,m1,n_old)/derivative2(z,k,xxi,m1,n_old);
		if (n_new<0) {
			n_new = 0.0;
		};
		if (n_new>1) {
			n_new = 1.0;
		};
		double error = abs((n_new-n_old)/n_new);
		if (error<mytol) {
			found = 1;
		} else {
			n_old = n_new;
		};

		iter++;
	};

	if (found==1) {
		return n_new;
	} else {
		return -99.99;
	};
};

// Eureka function check whether a tuple (STATE,SHOCK,SHADOW) can survive to next iteration
__host__ __device__
int eureka(double k, double z, double xi,
            double m1, int i_z, int i_xi,
            double* K, 
			double* EM1_low, double* EM1_high, double stuff
			) {

	// Declare Variables
	double d, mu, n, Y, c, kplus, w, lhs1, interp_low, interp_high;
	int i_kplus;

	// Case 1: Not Binding
	mu = 0;
	n = newtonlabor(z,k,xi,m1,1);
	Y = stuff*pow(n,1-ttheta);
	c = (1-ttheta)*Y/(aalpha*n);
	kplus = (1-ddelta)*k + Y - c;
	w = aalpha*c;
	d = c - w*n;

	lhs1 = (1-mu*xi)/c;

	i_kplus = fit2grid(kplus,nk,K);
	interp_low = EM1_low[i_kplus+nk*(i_z+i_xi*nz)] + 
		(EM1_low[i_kplus+1+nk*(i_z+i_xi*nz)]-EM1_low[i_kplus+nk*(i_z+i_xi*nz)])/(K[i_kplus+1]-K[i_kplus])*(kplus-K[i_kplus]);
	interp_high = EM1_high[i_kplus+nk*(i_z+i_xi*nz)] + 
		(EM1_high[i_kplus+1+nk*(i_z+i_xi*nz)]-EM1_high[i_kplus+nk*(i_z+i_xi*nz)])/(K[i_kplus+1]-K[i_kplus])*(kplus-K[i_kplus]);
	if (
		(bbeta*interp_low <= lhs1) &&
		(lhs1 <=bbeta*interp_high) &&
		(c>0) && (xi*kplus>w*n) && (kplus>0) && (n>0) && (n<=1)
	   )
	{
		// printf("Found  Not Binding\n");
		return 1;
	}; 


	// Case 2: Binding
	n = newtonlabor(z,k,xi,m1,2); 
	c = (1-ddelta)/m1 + ttheta*stuff*pow(n,1-ttheta)/(m1*k); // correct
	w = aalpha*c;	// correct
	kplus = w*n/xi;	// correct
	d = c - w*n;	// correct
	mu = (1-ttheta)*stuff*pow(n,-ttheta)/w - 1;	// correct

	i_kplus = fit2grid(kplus,nk,K);
	interp_low = EM1_low[i_kplus+nk*(i_z+i_xi*nz)] + 
		(EM1_low[i_kplus+1+nk*(i_z+i_xi*nz)]-EM1_low[i_kplus+nk*(i_z+i_xi*nz)])/(K[i_kplus+1]-K[i_kplus])*(kplus-K[i_kplus]);
	interp_high = EM1_high[i_kplus+nk*(i_z+i_xi*nz)] + 
		(EM1_high[i_kplus+1+nk*(i_z+i_xi*nz)]-EM1_high[i_kplus+nk*(i_z+i_xi*nz)])/(K[i_kplus+1]-K[i_kplus])*(kplus-K[i_kplus]);
	lhs1 = (1-mu*(xi))/c;

	if (
		(bbeta*interp_low <= lhs1) &&
		(lhs1 <=bbeta*interp_high) &&
		(c>0) && (kplus>0) && (mu>=0) && (n>0) && (n<=1)
	   )
	{
		// printf("Found Binding \n");
		return 2;
	};

	return 0;
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
	// Steady-State calculation
	double zss = 1;
	double xxiss = xxibar;
	double nomin = xxiss+1-bbeta*(1-ddelta);
	double denom = 1-ttheta+bbeta*ttheta; 
	double kovern = pow(nomin/denom,1/(ttheta-1));
	double covern = pow(kovern,ttheta) - ddelta*kovern;
	double css = xxiss*kovern/aalpha;
	double nss = css/covern;
	double kss = kovern*nss;
	double wss = aalpha*css;
	double mmuss = (1-bbeta*(1-ddelta+ttheta*pow(kovern,ttheta-1)))/xxiss;
	double dss = css - wss*nss;
	double invss = ddelta*kss;
	double mkss =(1/css)*( 1-ddelta+ttheta*zss*pow(kovern,ttheta-1)  );
	double yss = zss*pow(kss,ttheta)*pow(nss,1-ttheta);
	double lhsfc = (xxiss*kss);
	double rhsfc = wss*nss;

	cout << setprecision(16) << "kss: " << kss << endl;
	cout << setprecision(16) << "invss: " << invss << endl;
	cout << setprecision(16) << "zss: " << zss << endl;
	cout << setprecision(16) << "xxiss: " <<xxibar << endl;
	cout << setprecision(16) << "mkss: " << mkss << endl;
	cout << setprecision(16) << "dss: " << dss << endl;
	cout << setprecision(16) << "css: " << css << endl;
	cout << setprecision(16) << "nss: " << nss << endl;
	cout << setprecision(16) << "wss: " << wss << endl;
	cout << setprecision(16) << "mmuss: " << mmuss << endl;
	cout << setprecision(16) << "lhsfc: " << lhsfc << endl;
	cout << setprecision(16) << "rhsfc: " << rhsfc << endl;

	// Select Device
	int num_devices;
	cudaGetDeviceCount(&num_devices);
	if (num_devices>1) {
		cudaSetDevice(0);
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
	chebyroots_cuda(nk,h_K_cheby_ptr);
	h_K = h_K_cheby;
	double* h_K_ptr = raw_pointer_cast(h_K.data());
	double minK = (1-kwidth)*kss;
	double maxK = (1+kwidth)*kss;
	cout << "minK: " << minK << endl;
	cout << "maxK: " << maxK << endl;
	fromchebydomain(minK, maxK, nk, h_K_ptr);

	// initialize shock grids and transition matrix
	double A [4];
	A[0] = rrhozz;
	A[1] = rrhoxxiz;
	A[2] = rrhozxxi;
	A[3] = rrhoxxixxi;
	double Ssigma [4];
	Ssigma[0] = ssigmaz*ssigmaz;
	Ssigma[1] = 0;
	Ssigma[2] = 0;
	Ssigma[3] = ssigmaxxi*ssigmaxxi;
	double* h_P_ptr = raw_pointer_cast(h_P.data());
	host_vector<double> h_shockgrids(nz+nxxi);
	double* h_shockgrids_ptr = raw_pointer_cast(h_shockgrids.data());
	tauchen_vec(2, nz, 3, A, Ssigma, h_shockgrids_ptr, h_P_ptr);
	for (int index=0; index<nz; index++) {
		h_Z[index] = zbar*exp(h_shockgrids[index+0*nz]);
		h_XXI[index] = xxibar*exp(h_shockgrids[index+1*nz]);
	};
	double* h_Z_cheby_ptr = raw_pointer_cast(h_Z_cheby.data());
	double* h_XXI_cheby_ptr = raw_pointer_cast(h_XXI_cheby.data());
	chebyroots_cuda(nz,h_Z_cheby_ptr);
	chebyroots_cuda(nxxi,h_XXI_cheby_ptr);
	printstuff< host_vector<double> >(h_XXI);

	// Read the initial guess
	ifstream fin_guess("./MATLAB/coeff_guess.txt");
	for (int i=0; i<(1+pk)*(1+pz)*(1+pxxi); i++) {
		double temp;
		fin_guess >> setprecision(18) >> temp;
		h_coeff[i] = temp;
	};
	fin_guess.close();
	printstuff< host_vector<double> > (h_coeff);

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

	// Device only vectors
	// device_vector<double> d_X(nk*nz*nxxi*pk*pz*pxxi);
	// device_vector<double> d_M(nk*nz*nxxi);

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
	d_projector = h_projector;

	


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

		// printstuff< host_vector<double> > (d_M);

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

		// Find EMs for low and high 
		// cublasDgemm(handle,
		// 		CUBLAS_OP_N,  
		// 		CUBLAS_OP_T,
		// 		nk*nb, nz*nxi, nz*nxi,
		// 		&alpha,
		// 		d_V1_low_ptr, 
		// 		nk*nb, 
		// 		d_P_ptr,
		// 		nz*nxi,
		// 		&beta,
		// 		d_EM1_low_ptr,
		// 		nk*nb);
		// cublasDgemm(handle,
		// 		CUBLAS_OP_N,  
		// 		CUBLAS_OP_T,
		// 		nk*nb, nz*nxi, nz*nxi,
		// 		&alpha,
		// 		d_V1_high_ptr,  
		// 		nk*nb,      
		// 		d_P_ptr,
		// 		nz*nxi,
		// 		&beta,
		// 		d_EM1_high_ptr,
		// 		nk*nb);

		// // Directly find the new Value function
		// thrust::for_each(
		// 	begin,
		// 	end,
		// 	RHS(d_K_ptr, d_Z_ptr, d_XI_ptr,
		// 		d_V1_low_ptr,
		// 		d_V1_high_ptr,
		// 		d_Vplus1_low_ptr,
		// 		d_Vplus1_high_ptr,
		// 		d_EM1_low_ptr,
		// 		d_EM1_high_ptr,
		// 		d_flag_ptr)
		// );

		// // Find error
		// double diff1 = transform_reduce(
		// 	make_zip_iterator(make_tuple(d_V1_low.begin(), d_Vplus1_low.begin(), d_V1_high.begin(),d_Vplus1_high.begin())),
		// 	make_zip_iterator(make_tuple(d_V1_low.end()  , d_Vplus1_low.end()  , d_V1_high.end()  ,d_Vplus1_high.end())),
		// 	myMinus(),
		// 	0.0,
		// 	maximum<double>()
		// 	);

		// // Find distance 
		// double dist1 = transform_reduce(
		// 	make_zip_iterator(make_tuple(d_Vplus1_low.begin(),d_Vplus1_high.begin())),
		// 	make_zip_iterator(make_tuple(d_Vplus1_low.end()  ,d_Vplus1_high.end())),
		// 	myDist(),
		// 	0.0,
		// 	maximum<double>()
		// 	);
		// diff = max(diff1,-99.0);
		// dist = max(dist1,-99.0);

		// cout << "diff is: "<< diff << endl;
		// cout << "dist is: "<< dist << endl;
		// cout << "Vplus1[1035,10,6] (the spike) range is " << d_Vplus1_low[1035+nk*10+nk*nz*6] << ", " << d_Vplus1_high[1035+nk*10+nk*nz*6] << endl;

		// // update correspondence
		// d_V1_low = d_Vplus1_low; d_V1_high = d_Vplus1_high;

		// cout << ++iter << endl;
		// cout << "=====================" << endl;
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


