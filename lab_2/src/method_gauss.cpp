#include "method_gauss.h"

#include "Matrix_public.h"

#define eps 1e-16

#include <iomanip>
#include <cmath>
#include <iostream>

int method_gauss(double**& A, size_t& n, double*& x) {

	int check = 0;

	check = straight_course(A, n);

	if (check == 0)
		reverse_course(A, n, x);
	else
		return 1;

	return 0;
}

//---subfunctions1---//

int straight_course(double**& A, size_t& n) {
	//printMatrix(A, n);
	major_element(A, n);
	//printMatrix(A, n);
	for (size_t k = 0; k < n; ++k) {

		if (fabs(A[k][k]) < eps)
		{
			cout << "\nMatrix is singular (non-invertible)\n";
			return 1;

		}	
		remove_MatrixColumnElements_UnderLine(k, A, n);

		//print_Sys_Matrix(A, n);
	}
	return 0;

}

int reverse_course(double**& A, size_t& n, double*& x) {


	for (size_t i = 0; i < n; ++i) {
		//normalize_Sys_MatrixLine(i, A, n);
		x[i] = A[i][n];
	}
	x[n - 1] /= A[n - 1][n - 1];
	for (size_t i = 0; i <= n - 2; ++i) {

		for (size_t j = n - 2 - i + 1; j < n; ++j) {
			x[n - 2 - i] -= x[j] * A[n - 2 - i][j];

		}
		x[n - 2 - i] /= A[n - 2 - i][n - 2 - i];
	}

	check_vector_zero(x, n);
	
	return 0;
}

//---subfunctions2---//

int major_element(double**& A, size_t& n) {

	size_t i_max = 0;

	for (size_t k = 0; k < n; ++k) {
		i_max = k;
		for (size_t i = k; i < n; ++i) {
			if (fabs(fabs(A[i_max][k]) - fabs(A[i][k])) > eps) {
				if (fabs(A[i_max][k]) < fabs(A[i][k])) {
					i_max = i;
				}

			}

		}
		swap(A[i_max], A[k]);
	}
	return 0;
}

int remove_MatrixColumnElements_UnderLine(size_t& k, double**& A, size_t& n) {

	double tmp1 = 0.0;

	for (size_t i = k + 1; i < n; ++i) {


		if (fabs(A[i][k]) > eps) {

			tmp1 = A[i][k]/A[k][k];

			for (size_t j = k; j < n + 1; ++j) {

				A[i][j] -= A[k][j] * tmp1;
				//printMatrix(A, n);
				//solving -0 problem
				if (fabs(A[i][j]) < eps) {
					A[i][j] *= (1 / A[i][j] > 0) ? (1) : (-1);
				}
			}
		}
		else {
			A[i][k] = 0.0;
			if (fabs(A[i][k]) < eps) {
				A[i][k] *= (1 / A[i][k] > 0) ? (1) : (-1);
			};
		}

	}

	return 0;
}