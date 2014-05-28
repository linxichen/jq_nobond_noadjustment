// #include <thrust/host_vector.h> // This include happens in linking stage
#include <iostream>
#include <iomanip>

void linspace(double, double, int, double*);
void mynormalcpp(double *, double, double, int, unsigned);
void myuniformcpp(double *, int, unsigned);
void myexponentialcpp(double *, int, unsigned);
void mytest(const int);
void fromchebydomain(double, double, int, double*);
void ind2sub(int, int*, int, int*);
void findprojector(double*, int, int, double*);

typedef void (*gridgen_fptr)(double, double, int, double*);
void tauchen_vec(int, int, int, double*, double*, double*, double*, gridgen_fptr);

