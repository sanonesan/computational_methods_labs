#include "QR_method.h"
#include <cmath>
#include "Matrix_public.h"
#include "method_gauss.h"
#include <iostream>
#include <cmath>

#define eps 1e-16

int QR_method_find_Q_and_R(
	double** &Q, 
	double** &R, 
	double** &A, size_t& n, double*& x) {

	//print_Sys_Matrix(A, n);

	Rb_sys_matrix(
		Q, 
		R, 
		A, n);

	for (size_t k = 0; k < n; ++k) 
		if (fabs(A[k][k]) < eps)
		{
			std::cout << "\nMatrix is singular (non-invertible)\n";
			return 1;
		}
	//print_Sys_Matrix(A, n);
	reverse_course(A, n, x);

	return 0;
}

//---subfunctions---//
int Rb_sys_matrix(
	double**& Q,
	double**& R, 
	double**& A, size_t& n) {

	double c = 0.0, s = 0.0;

	
	marix_make_E(Q, n);
	marix_make_E(R, n);

	double a = 0.0, b = 0.0;

	for (size_t i = 0; i < n - 1; ++i) {

		for (size_t j = i + 1; j < n; ++j) {

			coefs(i, j, c, s, A, n);

			for (size_t k = 0; k < n+1; ++k) {

				a = A[i][k];
				b = A[j][k];
				A[i][k] = (c * a + s * b);
				A[j][k] = (-s * a + c * b);

				if (k < n) {
					a = Q[i][k];
					b = Q[j][k];
					Q[i][k] = (c * a + s * b);
					Q[j][k] = (-s * a + c * b);
				}
				
				
			}				

		}

	}
	

	copy_Matrix_A_to_B(A, R, n);

	check_matrix_zero(A, n);
	check_matrix_zero(Q, n);
	check_matrix_zero(R, n);




	return 0;
}

int coefs(size_t& k, size_t& l, double& c, double& s, double**& A, size_t& n) {

	if (k < n && l < n) {

		double temp = sqrt(pow(A[k][k], 2) + pow(A[l][k], 2));

		c = A[k][k] / temp;

		s = A[l][k] / temp;

		if (fabs(c) < eps)
			c = (1 / c > 0) ? c : c*(-1);
		if (fabs(s) < eps)
			s = (1 / s > 0) ? s : s*(-1);

	}


	return 0;
}

