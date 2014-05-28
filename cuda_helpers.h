#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

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

// A function template to display vectors, C array style
template <class T>
__host__ __device__
void display_vec(T vec, int size) {
	for (int i = 0; i < size; i++) {
		std::printf("The %ith element, @[%i] = %f\n", i+1, i, vec[i]);
	};
};

// A function template to display vectors, std::vector style
template <class T>
__host__ __device__
void display_vec(T vec) {
	int size = vec.size();
	for (int i = 0; i < size; i++) {
		std::printf("The %ith element, @[%i] = %f\n", i+1, i, vec[i]);
	};
};

// A function template to save vectors to file, C array style
template <class T>
void save_vec(T vec, int size, std::string filename ) {
	std::ofstream fileout(filename.c_str(), std::ofstream::trunc);
	for (int i = 0; i < size; i++) {
		fileout << vec[i] << '\n';
	};
	fileout.close();

};

// A function template to save vectors to file, std::vector style
template <class T>
void save_vec(T vec, std::string filename ) {
	int size = vec.size();
	std::ofstream fileout(filename.c_str(), std::ofstream::trunc);
	for (int i = 0; i < size; i++) {
		fileout << vec[i] << '\n';
	};
	fileout.close();
};

// Newton's Method with bracketing, i.e. we know on two points the function differs in sign.
// Codes from Numerical Recipes 3rd. Ed.
// BEWARE: The stopping criteria is not right yet.
template <class T>
__host__ __device__
double newton_bracket(T func, const double x1, const double x2, double x0) {
// Purpose: Tries to find a root for function named func. Its first derivative is given by func.prime().
//			It is assumed that func(x1) and func(x2) are different in sign so a root exists within. x0 is the guess.
	const int maxiter = 100;
	const double tol = 1e-5;
	// Checking the bounds: they need to make sense. Or sometimes the bounds are solutions.
	double f1 = func(x1);
	double f2 = func(x2);
	if (f1*f2>0) return -5179394.1; // The different sign assumption violated!
	if (f1 == 0) return x1;
	if (f2 == 0) return x2;

	// Orient the search so that f(xl) < 0
	double xl, xh;
	if (f1 < 0.0) {
		xl = x1;
		xh = x2; 
	} else {
		xh = x1;
		xl = x2;
	};

	// Initialize guess and other things
	double rts = x0;
	double dxold = abs(x2-x1);
	double dx = dxold;
	double f = func(rts);
	double df = func.prime(rts);

	for (int iter = 0; iter < maxiter; iter++) { 
		if ( 
			( ((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0 )   ||	// Bisect if Newton step out of range
			( abs(2.0*f) > abs(dxold*df)  ) // ... or step not decreasing fast enough
		   )
		{
			dxold = dx;
			dx = 0.5*(xh-xl);
			rts += dxold; // undo the newton step
			rts = xl + dx;
			if (xl == rts) return rts;
		} else {
			// If newton step is okay 
			dxold = dx;
			dx = f/df;
			double temp = rts;
			rts -= dx;
			if (temp==rts) return rts;
		};

		// Check for convergence
		if ( abs(dx)/(1+abs(rts+dx)) < tol ) return rts;

		// Compute new f and df for next iteration
		f = func(rts);
		df = func.prime(rts);

		// Maintain the bracket
		if (f < 0.0) {
			xl = rts;
		} else {
			xh = rts;
		};
	};

	return -51709394.2;
};

// "Raw" Newton's Method
// Codes from Numerical Recipes 3rd. Ed.
template <class T>
__host__ __device__
double newton(T func, const double x1, const double x2, double x0) {
// Purpose: Tries to find a root for function named func. Its first derivative is given by func.prime().
//			func is only defined on [x1,x2] We "pull back" when outside. x0 is the guess.
	const int maxiter = 20;
	const double tol = 1e-5;
	// Initialize guess and other things
	double x_old = x0;
	double x = x0;
	for (int iter = 0; iter < maxiter; iter++) { 
		double f1 = func(x1);
		double f2 = func(x2);
		if (f1==0) return x1;
		if (f2==0) return x2;
		x = x_old - func(x)/func.prime(x);

		// Pull back if outside of support
		if (x<x1) {
			x = x1;
		};
		if (x>x2) {
			x = x2;
		};

		// Check for convergence
		if ( (abs(x-x_old)/(1+abs(x_old))<tol) && (x>x1) && (x<x2) ) {
		   	return x;
		} else {
			x_old = x;
		};
	};

	return -51709394.2;
};
