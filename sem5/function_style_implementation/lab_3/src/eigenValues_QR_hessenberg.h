
#pragma once 

#include <iostream>
#include <cmath>

#define eps 1e-16

using namespace std;

double* hessenberg(double**& A, size_t& n);

int coefs_hessen(size_t& k, size_t& l, double& alpha, double& beta, double**& A, size_t& n);