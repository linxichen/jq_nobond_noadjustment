/* This is a CUDA implementation of Jermann Quadrini 2013 AER  
 * We simulate a equilibrium that collateral constraint is always binding to check the accuracy of 
 * their linearization approach. Hopeully we can ind something that they missed. A main suspect 
 * is the asymmetry of policy functions.
 */

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
// #include <cuda_runtime.h>
#include <cublas_v2.h>
// #include <helper_functions.h>

// Includes, C++ codes
#include "cppcode.h"

using namespace std;
using namespace thrust;

// Accuracy controls
// #define nk 256
// #define nb 1 
// #define nz 3
// #define nxi 3
// #define nm1 256 
// #define tol 1e-10
// #define maxiter 1
// #define kwidth 2.0 
// #define bwidth 1.15 
// #define mkwidth 15.0 // has to be a double 

// // Parameters from the AER paper
// #define aalpha 1.8834
// #define bbeta 0.9825
// #define ddelta 0.025
// #define ttheta 0.36
// #define kkappa 0.1460
// #define ttau 0.3500
// #define xxibar 0.1634
// #define zbar 1.0
// #define dbar 2.56514

// Define an class that contains parameters
struct para_struct {
	// Accuracy controls
	int nk ;
	int nb ;
	int nz ;
	int nxxi;
	int nm1;
	int tol;
	int maxiter;
	double kwidth ;
	double bwidth ;
	double mkwidth;

	// Model parameters
	int aalpha;
	int bbeta ;
	int ddelta;
	int ttheta;
	int kkappa;
	int ttau  ;
	int xxibar;
	int zbar  ;
	int dbar  ;
};

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

/*********************************************************************************
// The objective function for case 2
__host__ __device__
double objective2(double z,double k,double xxi,double m1,double n,para_struct para) {
	return (1-para.ddelta)*k-(1-para.ddelta)/m1-(1-para.ddelta)*para.aalpha*n/(xxi*m1)
		+ z*pow(k,para.ttheta)*pow(n,1-para.ttheta)*(1-para.ttheta/(m1*k)-para.ttheta*para.aalpha*n/(m1*xxi*k));
};

// The  derivative for case 2
__host__ __device__
double derivative2(double z,double k,double xxi,double m1,double n,para_struct para) {
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
********************************************************************************/

// Define the type function ptr that takes four double and returns one doulbe
typedef double (* func_ptr)(double, double ,double, double, double);

// Apply Newton's Method to find labor hours. Some ideas are from the book "Numerical Recipes"
__host__ __device__
double newtonlabor(double z,double k,double xxi,double m1, func_ptr obj_func, func_ptr drv_func) {
	double n_old = 0.3;
	double n_new = 0.3;
	double mytol = 0.0000001;
	int mymaxiter = 30;
	int found = 0; // 1 means we found the solution

	int iter = 0;
	while ( (found != 1) && (iter <= mymaxiter) ) {
		n_new = n_old - obj_func(z,k,xxi,m1,n_old)/drv_func(z,k,xxi,m1,n_old);
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

	if ((found==1) && (iter<mymaxiter)) {
		return n_new;
	} else {
		return -99.99;
	};
};


// Eureka function check whether a tuple (STATE,SHOCK,SHADOW) can survive to next iteration
// now M4 is for stock price !!!!!!!
__host__ __device__
bool eureka(double k, double z, double xi,
            double m1, int i_z, int i_xi,
            double* K, 
			double* EM1_low, double* EM1_high, 
			double MPK, double Y,
			para_struct para) {

	// A test call of newtonlabor
	func_ptr obj1_ptr = & objective1;
	func_ptr drv1_ptr = & derivative1;
	double test = newtonlabor(z,k,xi,m1,obj1_ptr,drv1_ptr);

	// Case 1: Binding
	double kplus = Y/xi;
	double c = Y+(1-para.ddelta)*k-kplus;
	double mu = 1-(m1*c-1+para.ddelta)/MPK;	
	double w = (1-mu)*(1-para.ttheta)*Y;
	double d = c-w;

	double R = 1;
	int i_kplus = fit2grid(kplus,para.nk,K);
	// int i_bplus = fit2grid(bplus,nb,B);
	double lhs1 = (1-xi*mu)/c;
	
	if (
		(para.bbeta*EM1_low[i_kplus+para.nk*(i_z+i_xi*para.nz)] <= lhs1) &&
		(lhs1 <=para.bbeta*EM1_high[i_kplus+para.nk*(i_z+i_xi*para.nz)]) &&
		(c>0) && (mu>=0) && (w>=0) //  &&  
		// (kplus <= K[nk-1]) && (kplus >= K[0])
	   )
	{
		return true;
	};

	// Case 2: Not Binding
	mu = 0;
	c = (1-para.ddelta+MPK)/m1;	
	kplus = (1-para.ddelta)*k + Y - c;
	w = (1-para.ttheta)*Y;
	d = c - w;
	lhs1 = (1-xi*mu)/c;
	i_kplus = fit2grid(kplus,para.nk,K);
	if (
		(para.bbeta*EM1_low[i_kplus+para.nk*(i_z+i_xi*para.nz)] <= lhs1) &&
		(lhs1 <=para.bbeta*EM1_high[i_kplus+para.nk*(i_z+i_xi*para.nz)]) &&
		(c>0) && (w>=0) && (xi*kplus>Y) // &&
		// (kplus <= K[nk-1]) && (kplus >= K[0])
	   )
	{
		return true;
	} else {
		return false;
	};
};

// This functor yields RHS for each (k', k, z). Follwing examples in Thrust doc
struct RHS 
{
	// Data member
	double *K, *Z, *XI;
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
	RHS(double* K_ptr, double* Z_ptr, double* XI_ptr,
	double* V1_low_ptr,
	double* V1_high_ptr,
	double* Vplus1_low_ptr,
	double* Vplus1_high_ptr,
	double* EM1_low_ptr,
	double* EM1_high_ptr,
	double* flag_ptr,
	para_struct _para)
	{
		K = K_ptr; Z = Z_ptr; XI = XI_ptr;
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
		int i_xi = index/(para.nk*para.nz);
		int i_z  = (index-i_xi*para.nk*para.nz)/(para.nk);
		int i_k = index - i_xi*para.nk*para.nz - i_z*para.nk;
		double k =K[i_k]; double z=Z[i_z]; double xi=XI[i_xi];
		double Y = z*pow(k,para.ttheta); double MPK = para.ttheta*z*pow(k,para.ttheta-1);

		// Find the "box" or "hypercube" that described m's range. Fancy word.
		double m1min = V1_low[index]; double m1max = V1_high[index];
		double step1 = (m1max-m1min)/double(para.nm1-1);
		double tempflag = 0.0;

		// Brute force searching
		int i_m1_min=0;
		int i_m1_max=para.nm1-1;
		int resume_ind=0; //
		// Initial search to have something to compare later
		for (int m_index = 0; m_index < para.nm1; m_index++) {
			int i_m1 = m_index;
			double m1=m1min+double(i_m1)*step1;
			if (eureka(k,z,xi,m1,i_z,i_xi,K,
				EM1_low,EM1_high,MPK,Y))
			{
				tempflag++;
				i_m1_max = i_m1;
				i_m1_min = i_m1;
				resume_ind = m_index + 1;
				break; // break and only break one layer of for loop
			};
		};

		// Update min and max, run to the end
		for (int m_index = resume_ind; m_index < para.nm1; m_index++) {
			int i_m1 = m_index;
			double m1=m1min+double(i_m1)*step1;
			if (eureka(k,z,xi,m1,i_z,i_xi,K,
				EM1_low,EM1_high,MPK,Y))
			{
				tempflag++;
				if (i_m1 > i_m1_max) i_m1_max = i_m1;
				if (i_m1 < i_m1_min) i_m1_min = i_m1;
			};
		};

		// Update Vs
		flag[index] = double(tempflag)/double(para.nm1);
		if (tempflag == 0) {
			Vplus1_high[index] = 1*(m1max+0.0*step1);
			Vplus1_low[index] = 1*(m1min-0.0*step1);
		} else {
			Vplus1_high[index] = m1min+i_m1_max*step1;
			Vplus1_low[index] = m1min+i_m1_min*step1;
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

// Main
int main()
{
	// Set Parameters
	para_struct para;

	// Set Accuracy Parameters
	para.nk = 256;
	para.nb = 1 ;
	para.nz = 3;
	para.nxi = 3;
	para.nm1 = 256 ;
	para.tol = 1e-10;
	para.maxiter = 1;
	para.kwidth = 2.0 ;
	para.bwidth = 1.15 ;
	para.mkwidth = 15.0 ; 

	para.aalpha = 1.8834;
	para.bbeta = 0.9825;
	para.ddelta = 0.025;
	para.ttheta = 0.36;
	para.kkappa = 0.1460;
	para.ttau = 0.3500;
	para.xxibar = 0.1634;
	para.zbar = 1.0;
	para.dbar = 2.56514;

	// Steady-State calculation
	double zss = 1;
	double nss = 1;
		
	double mmuss = 0;
	double kss = pow((1/para.bbeta-1+para.ddelta)/para.ttheta,1/(para.ttheta-1));
	double css = pow(kss,para.ttheta)+(1-para.ddelta)*kss-kss;
	double wss = (1-para.ttheta)*pow(kss,para.ttheta);
	double dss = css - wss;
	double pss = dss/(1-para.bbeta);
	double lhsfc = para.xxibar*kss;
	double rhsfc = pow(kss,para.ttheta);
	double mkss = (1-para.ddelta+para.ttheta*pow(kss,para.ttheta-1))/css;
	double msss = pss/css;
	
	cout << setprecision(16) << "kss: " << kss << endl;
	cout << setprecision(16) << "zss: " << para.zbar << endl;
	cout << setprecision(16) << "xxiss: " <<para.xxibar << endl;
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
	host_vector<double> h_K(para.nk); 
	host_vector<double> h_Z(para.nz);
	host_vector<double> h_XI(para.nxxi);

	host_vector<double> h_V1_low(para.nk*para.nz*para.nxxi, 1/para.mkwidth*mkss);
	host_vector<double> h_V1_high(para.nk*para.nz*para.nxi,para.mkwidth*mkss);
	host_vector<double> h_Vplus1_low(para.nk*para.nz*para.nxxi,1/para.mkwidth*mkss);
	host_vector<double> h_Vplus1_high(para.nk*para.nz*para.nxxi,para.mkwidth*mkss);

	host_vector<double> h_EM1_low(para.nk*para.nz*para.nxxi,0.0);
	host_vector<double> h_EM1_high(para.nk*para.nz*para.nxxi,0.0);

	host_vector<double> h_P(para.nz*para.nxxi*para.nz*para.nxxi, 1.0/double(para.nz*para.nxxi)); // No Tauchen yet
	host_vector<double> h_flag(para.nk*para.nz*para.nxxi, 0); 
	
	// initialize capital grid
	double minK = 1/para.kwidth*kss;
	double maxK = para.kwidth*kss;
	double step = (maxK-minK)/double(para.nk-1);
	for (int i_k = 0; i_k < para.nk; i_k++) {
		h_K[i_k] = minK + step*double(i_k);
	};

	// initialize TFP shock grid
	double minZ = 0.8*para.zbar;
	double maxZ = 1.2*para.zbar;
	step = (nz>1)?(maxZ - minZ)/double(nz-1): 0;
	for (int i_z = 0; i_z < nz; i_z++) {
		h_Z[i_z] = minZ + step*double(i_z);
	};

	// initialize financial shock grid
	double minXI = 0.8*para.xxibar;
	double maxXI = 1.2*para.xxibar;
	step = (nxi>1)? (maxXI - minXI)/double(para.nxxi-1): 0;
	for (int i_xi = 0; i_xi < para.nxxi; i_xi++) {
		h_XI[i_xi] = minXI + step*double(i_xi);
	};

	// Copy to the device
	device_vector<double> d_K = h_K;
	device_vector<double> d_Z = h_Z;
	device_vector<double> d_XI = h_XI;

	device_vector<double> d_V1_low = h_V1_low;
	device_vector<double> d_V1_high = h_V1_high;

	device_vector<double> d_Vplus1_low = h_Vplus1_low;
	device_vector<double> d_Vplus1_high = h_Vplus1_high;

	device_vector<double> d_EM1_low = h_EM1_low;
	device_vector<double> d_EM1_high = h_EM1_high;

	device_vector<double> d_P = h_P;
	device_vector<double> d_flag = h_flag;
	// device_vector<double> d_keys_out(nk*nz);

	// Obtain device pointers to be used by cuBLAS
	double* d_K_ptr = raw_pointer_cast(d_K.data());
	double* d_Z_ptr = raw_pointer_cast(d_Z.data());
	double* d_XI_ptr = raw_pointer_cast(d_XI.data());

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
	counting_iterator<int> end(nk*nz*nxi);

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
	while ((diff>tol)&&(iter<maxiter)){
		// Find EMs for low and high 
		cublasDgemm(handle,
				CUBLAS_OP_N,  
				CUBLAS_OP_T,
				nk*nb, nz*nxi, nz*nxi,
				&alpha,
				d_V1_low_ptr, 
				nk*nb, 
				d_P_ptr,
				nz*nxi,
				&beta,
				d_EM1_low_ptr,
				nk*nb);
		cublasDgemm(handle,
				CUBLAS_OP_N,  
				CUBLAS_OP_T,
				nk*nb, nz*nxi, nz*nxi,
				&alpha,
				d_V1_high_ptr,  
				nk*nb,      
				d_P_ptr,
				nz*nxi,
				&beta,
				d_EM1_high_ptr,
				nk*nb);

		// cublasDgemm(handle,
		// 		CUBLAS_OP_N,  
		// 		CUBLAS_OP_T,
		// 		nk*nb, nz*nxi, nz*nxi,
		// 		&alpha,
		// 		d_V4_low_ptr, 
		// 		nk*nb, 
		// 		d_P_ptr,
		// 		nz*nxi,
		// 		&beta,
		// 		d_EM4_low_ptr,
		// 		nk*nb);
		// cublasDgemm(handle,
		// 		CUBLAS_OP_N,  
		// 		CUBLAS_OP_T,
		// 		nk*nb, nz*nxi, nz*nxi,
		// 		&alpha,
		// 		d_V4_high_ptr,  
		// 		nk*nb,      
		// 		d_P_ptr,
		// 		nz*nxi,
		// 		&beta,
		// 		d_EM4_high_ptr,
		// 		nk*nb);

		// Directly find the new Value function
		thrust::for_each(
			begin,
			end,
			RHS(d_K_ptr, d_Z_ptr, d_XI_ptr,
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
	
	ofstream fout_V1_low("V1_low.csv", ios::trunc); ofstream fout_V1_high("V1_high.csv", ios::trunc);
	ofstream fout_kopt("koptcuda.csv", ios::trunc); ofstream fout_copt("coptcuda.csv", ios::trunc);
	ofstream fout_R("Rcuda.csv", ios::trunc);
	ofstream fout_d("dcuda.csv", ios::trunc); ofstream fout_n("ncuda.csv", ios::trunc);
	ofstream fout_Kgrid("Kgridcuda.csv", ios::trunc);
	ofstream fout_Zgrid("Zgridcuda.csv", ios::trunc); ofstream fout_XIgrid("XIgridcuda.csv", ios::trunc);
	ofstream fout_mmu("mmucuda.csv", ios::trunc); ofstream fout_P("Pcuda.csv", ios::trunc);
	ofstream fout_flag("flagcuda.csv", ios::trunc);
	
	for (int index=0; index<nk*nb*nz*nxi; index++) {
		fout_V1_low << h_V1_low[index] << '\n';
		fout_V1_high << h_V1_high[index] << '\n';
		int i_xi = index/(nk*nb*nz);
		int i_z  = (index-i_xi*nk*nb*nz)/(nk*nb);
		// int i_b  = (index-i_xi*nk*nb*nz-i_z*nk*nb)/nk;
		int i_k = index - i_xi*nk*nz - i_z*nk ;
		double m1 = (h_V1_low[index]+h_V1_low[index])/2;
		double k =h_K[i_k];
		// double b =h_B[i_b];
		double z=h_Z[i_z]; double xi=h_XI[i_xi];
		double Y = z*pow(k,ttheta);
		double MPK = ttheta*pow(k,ttheta-1);

		// Case 1: Binding
		double kplus = Y/xxibar;
		double c = Y+(1-ddelta)*k-kplus;
		double mu = 1-(m1*c-1+ddelta)/MPK;	
		double w = (1-mu)*(1-ttheta)*Y;
		double d = c-w;

		// int i_bplus = fit2grid(bplus,nb,B);
	
		if ( (mu>=0) )
		{
			fout_kopt << kplus << '\n';
			fout_copt << c << '\n';
			fout_R << 1 << '\n';
			fout_d << d << '\n';
			fout_n << 1 << '\n';
			fout_mmu << mu << '\n';
		} else {
			mu = 0;
			c = (1-ddelta+MPK)/m1;	
			kplus = (1-ddelta)*k + Y - c;
			w = (1-ttheta)*Y;
			d = c - w;
			fout_kopt << kplus << '\n';
			fout_copt << c << '\n';
			fout_R << 1 << '\n';
			fout_d << d << '\n';
			fout_n << 1 << '\n';
			fout_mmu << mu << '\n';
		};

		fout_flag << h_flag[index] << '\n';
	};
	
	for (int i=0; i<nk; i++) {
		fout_Kgrid << h_K[i] << '\n';
	};
	// for (int i=0; i<nb; i++) {
	// 	fout_Bgrid << h_B[i] << '\n';
	// };
	for (int i=0; i<nz; i++) {
		fout_Zgrid << h_Z[i] << '\n';
	};
	for (int i=0; i<nxi; i++) {
		fout_XIgrid << h_XI[i] << '\n';
	};
	for (int i=0; i<nz*nxi*nz*nxi; i++) {
		fout_P << h_P[i] << '\n';
	};

	fout_V1_low.close(); fout_V1_high.close();
	fout_Kgrid.close();
	fout_Zgrid.close(); fout_XIgrid.close();
	fout_kopt.close(); fout_copt.close();
	fout_R.close(); fout_d.close();
	fout_n.close(); fout_mmu.close(); fout_P.close();
	fout_flag.close();
	return 0;
}


