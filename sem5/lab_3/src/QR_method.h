#pragma once
#include <cstddef>
#define eps 1e-16


using namespace std;

int QR_method_find_Q_and_R(
	double** &Q, 
	double** &R, 
	double** &A, size_t& n, double*& x);

//---subfunctions---//
int Rb_sys_matrix(
	double** &Q, 
	double** &R, 
	double** &A, size_t& n);

int coefs(size_t &k, size_t &l, double &c, double &s, double** &A, size_t& n);

int QR_matrix(
	double**& Q,
	double**& R, 
	double**& A, size_t& n);






