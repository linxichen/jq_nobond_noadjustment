// #include <thrust/host_vector.h> // This include happens in linking stage
#include <iostream>
#include <iomanip>

void mynormalcpp(double *, double, double, int, unsigned);
void myuniformcpp(double *, int, unsigned);
void myexponentialcpp(double *, int, unsigned);
void mytest(const int);
void tauchen_vec(int, int, int, double*, double*, double*, double*);
void fromchebydomain(double, double, int, double*);
void ind2sub(int, int*, int, int*);
void findprojector(double*, int, int, double*);
