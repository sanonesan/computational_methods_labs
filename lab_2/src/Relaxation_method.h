#pragma once 
#include <cstddef>
#define eps 1e-16

double* Relaxation_method(double**& A, double*& x0, double& w, size_t n);

double* Relaxation_method(double*& a, double*& b, double*& c, double*& d, double*& x0, double& w, size_t n);
