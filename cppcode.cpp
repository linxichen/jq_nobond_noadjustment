#include <random>
#include <iostream>
#include <thrust/host_vector.h>
#include "armadillo"

using namespace arma;
using namespace std;

void mynormalcpp(double* output, double mmu, double ssigma, int n, unsigned seed) {
	
	// Specify engine and distribution
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(mmu,ssigma);

	// Actually generate random numbers
	for (int i=0; i<n; i++) {	
		output[i] = distribution(generator);
	};
	
};

void myuniformcpp(double* output, int n, unsigned seed)
{
	std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> d(0,1);
    
	for (int i=0; i<n; i++) {
		output[i] = d(generator);
	};
};

void myexponentialcpp(double* output, int n, unsigned seed)
{
	std::default_random_engine generator(seed);
    std::exponential_distribution<double> d(1.0);
    
	for (int i=0; i<n; i++) {
		output[i] = d(generator);
	};
};

int chebyroots_cpp(const int p, double* roots) {
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

void fromchebydomain(double lb, double ub, int p, double* ptr) {
	for (int i=0; i<p; i++) {
		ptr[i] = lb + (ptr[i]+1)/(2)*(ub-lb);	
	};
};

double mynormalcdf(double mean, double variance, double x) {
	return 0.5*( 1+erf( (x-mean)/sqrt(2*variance) ) );
};

// This function converts index to subscripts like ind2sub in MATLAB
void ind2sub(int length_size, int* siz_vec, int index, int* subs) {
// Purpose:		Converts index to subscripts. i -> [i_1, i_2, ..., i_n]
//
// Input:		length_size = # of coordinates, i.e. how many subscripts you are getting
// 				siz_vec = vector that stores the largest coordinate value for each subscripts. Or the dimensions of matrices
// 				index = the scalar index
//
// Ouput:		subs = the vector stores subscripts
	int cumdim [length_size];
	cumdim[0] = 1;
	for (int i=1; i<length_size; i++) {
		cumdim[i] = cumdim[i-1]*siz_vec[i-1];
	};
	int done = 0;
	for (int i=length_size-1; i>=0; i--) {
		subs[i] = (index - done)/cumdim[i];
		done += subs[i]*cumdim[i];
	};
};

// This function follows notation from Tauchen 1986 the vector case using Armadillo
void tauchen_vec(int M, int N, int m, double* A_ptr, double* Ssigma_e_ptr, double* Z_ptr, double* P_ptr) {
// Purpose:		Discretize y_t = A*y_t-1 + epsilon_t, where var(epsilon_t)=Ssigma_e
// 				y_t is an M by 1 vector
//
// Input:		M = # of shocks
// 				N = # of grid points for each shock
// 				m = # of s.d. away from mean should we use as bounds 
// 				A = Autocorrelation matrix
// 				Ssigma_e = variance-convariance matrix
//
// Output:		Z = Vectorized N by M matrix where each column stores grids
// 				P = Vectorized transition matrix (N^M by N^M)
	
	// 1. Read in matrices
	mat A(M,M);
	mat Ssigma_e(M,M);
	for (int i=0; i < M*M; i++) {
		int i_c = i/M;
		int i_r = i - i_c*M;
		A(i_r,i_c) = A_ptr[i];
		Ssigma_e(i_r,i_c) = Ssigma_e_ptr[i];
	};

	// 2. Find Variance Matrix for y
	vec Ssigma_e_vec = vectorise(Ssigma_e);
	vec Ssigma_y_vec = inv( eye(M*M,M*M) - kron(A,A) ) * Ssigma_e_vec;
	mat Ssigma_y = reshape(Ssigma_y_vec,M,M);
	A.print("Autocorr Matrix:");
	Ssigma_y.print("Var-Cov Matrix of Observables:");

	// 3. Create Grids. Assuming same # of grid points for each shock
	mat grids(M,N); // collects grid points here
	 
		// Option 1: Just Linspace
		// for (int i_shock=0; i_shock<M; i_shock++) {
		// 	double minshock = -m*sqrt(Ssigma_y(i_shock,i_shock));
		// 	double maxshock = +m*sqrt(Ssigma_y(i_shock,i_shock));
		// 	grids.row(i_shock) = trans(linspace(minshock,maxshock,N));
		// };
		// grids.print("Grids:");

		// Option 2: Chebyshev Nodes
		for (int i_shock=0; i_shock<M; i_shock++) {
			double minshock = -m*sqrt(Ssigma_y(i_shock,i_shock));
			double maxshock = +m*sqrt(Ssigma_y(i_shock,i_shock));
			cout << "minshock is: " << minshock << endl;
			cout << "maxshock is: " << maxshock << endl;
			vec temp_roots(N);
			chebyroots_cpp(N, temp_roots.memptr());
			fromchebydomain(minshock, maxshock, N, temp_roots.memptr());
			grids.row(i_shock) = trans(temp_roots);
		};
		grids.print("Grids:");

	// 4. Compute Transition Matrix
	mat P(pow(N,M),pow(N,M));
	cube h(M,pow(N,M),N);
	Col<int> sizes(M); sizes.fill(N);
	Col<int> subs(M);
	vec lag_y(M);
	for (int j=0; j<pow(N,M); j++) {
		// First find the conditional mean
		ind2sub(M, sizes.memptr(), j, subs.memptr());
		for (int i_shock=0; i_shock<M; i_shock++) lag_y(i_shock) = grids(i_shock,subs(i_shock));
		vec mmu = A*lag_y;

		// Fill the h cube
		double w_right, w_left;
		for (int i=0; i<M; i++) {
			for (int l=0; l<N; l++) {
				if (l==0) {
					w_right = grids(i,l+1) - grids(i,l);
					h(i,j,l) = mynormalcdf(0,Ssigma_e(i,i),grids(i,l)-mmu(i)+w_right/2);
				}
			   	else if (l==N-1) {
					w_left = grids(i,l) - grids(i,l-1);
					h(i,j,l) = 1 - mynormalcdf(0,Ssigma_e(i,i),grids(i,l)-mmu(i)-w_left/2);
				}
			   	else {
					w_left = grids(i,l) - grids(i,l-1);
					w_right = grids(i,l+1) - grids(i,l);
					h(i,j,l) = mynormalcdf(0,Ssigma_e(i,i),grids(i,l)-mmu(i)+w_right/2)
							 - mynormalcdf(0,Ssigma_e(i,i),grids(i,l)-mmu(i)-w_left/2);
				};
			};
		};
	};
	for (int j=0; j<pow(N,M); j++) {
		for (int k=0; k<pow(N,M); k++) {
			ind2sub(M, sizes.memptr(), k, subs.memptr());
			double cumprod = 1;
			for (int i=0; i<M; i++) {
				cumprod *= h(i,j,subs(i));
			};
			P(j,k) = cumprod;
		};
	};
	vec rowsum = sum(P,1);
	
	// 5. Check Accuracy
	arma_rng::set_seed(51709394);
	int T = 1e3;
	mat eps(M,T);
	for (int i_shock=0; i_shock<M; i_shock++) { 
		eps(i_shock,span::all) = trans(sqrt(Ssigma_e(i_shock,i_shock))*randn<vec>(T));
	};
	mat sim_y = zeros<mat>(M,T);
	for (int t=0; t<T-1; t++) {
		sim_y(span::all,t+1) = A*sim_y(span::all,t) + eps(span::all,t+1);
	};
	mat Z = sim_y(span::all,span(0,T-2));
	mat Y = sim_y(span::all,span(1,T-1));
	mat B(M,M);
	B = Y*trans(Z)*inv(Z*trans(Z)); //Multivariate Regression, formulas from Wikipedia
	B.print("Discrete Autocorr Matrix:");
	mat C = cov(trans(sim_y));
	C.print("Discrete Var-Cov Matrix of Observables: ");

	// 6. Output Results
	mat Z_out = trans(grids);
	double *Z_out_ptr = Z_out.memptr();
	for (int index=0; index<N*M; index++) {
		Z_ptr[index] = Z_out_ptr[index];
	};
	double* P_out_ptr = P.memptr();
	for (int index=0; index<pow(N,M)*pow(N,M); index++) {
		P_ptr[index] = P_out_ptr[index]; 
	};
};

void findprojector(double* X_ptr, int nrows, int ncols, double * P_ptr) {
	mat X(X_ptr, nrows, ncols, false, false);
	mat XprimeX = diagmat(trans(X)*X);
	// XprimeX.print("X'X=");
	mat P = inv(XprimeX)*trans(X);
	// P.print("Projector=");
	
	for (int i=0; i<nrows*ncols; i++) {
		P_ptr[i] = P.memptr()[i];
	};

};
