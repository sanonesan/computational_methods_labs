#include "QR_method.h"
#include <cmath>
#include "Matrix_public.h"
#include "method_gauss.h"
#include <iostream>




using namespace std;

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

	

	matrix_make_E(Q, n);
	matrix_make_E(R, n);

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



int QR_matrix(
	double**& Q,
	double**& R, 
	double**& A, size_t& n) {

	double c = 0.0, s = 0.0;


	double** A_save = new double* [n];

	for (size_t i = 0; i < n; ++i){
		A_save[i] = new double [n];
	}

	copy_Matrix_A_to_B(A, A_save, n);
	
	
	matrix_make_E(Q, n);
	matrix_make_E(R, n);

	double a = 0.0, b = 0.0;
	int iter = 0;
	for (size_t i = 0; i < n - 1; ++i) {

		for (size_t j = i + 1; j < n; ++j) {

			coefs(i, j, c, s, A, n);

			for (size_t k = 0; k < n; ++k) {

				a = A[i][k];
				b = A[j][k];

				if (a == 0 && b != 0){
					A[i][k] = (s * b);
					A[j][k] = (c * b);
					iter += 2;
				}
				if (a != 0 && b == 0){
					A[i][k] = (c * a);
					A[j][k] = (-s * a );
					iter += 2;

				}

				if (a == 0 && b == 0){
					A[i][k] = 0;
					A[j][k] = 0;
				}

				if (a != 0 && b != 0){
					
					A[i][k] = (c * a + s * b);
					A[j][k] = (-s * a + c * b);
					iter += 4;


				}



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

	copy_Matrix_A_to_B(A_save, A, n);

	check_matrix_zero(Q, n);

	transpose_matrix(Q, n);
	check_matrix_zero(R, n);

	for(size_t j = 0; j < n; ++j){
		delete [] A_save[j];
	}
	delete [] A_save;


	return iter;
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

