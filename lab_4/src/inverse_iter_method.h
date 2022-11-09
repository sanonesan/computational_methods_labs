#pragma once 
#include <cstddef>
#define eps 1e-15

double* inverse_iter_method(double**& A, double*& x0, double lambda, size_t n);


double* inverse_iter_method_Rayleigh(double**& A, double*& x0, double lambda, size_t n);
