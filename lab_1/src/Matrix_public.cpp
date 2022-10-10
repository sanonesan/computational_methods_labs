#include "Matrix_public.h"
#include "method_gauss.h"
#include "QR_method.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "time.h"
#include "stdlib.h"

using namespace std;


//--------------------//
double** readMatrixFromFile(string path, size_t& n) {

	ifstream file;
	file.open(path);

	if (file.is_open()) {

		file >> n;
		file.get();

		double** A = new double* [n];
		for (size_t j = 0; j < n; ++j)
			A[j] = new double[n + 1];

		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < n + 1; ++j) {
				file >> A[i][j];
			}

		file.close();

		return A;
	}
	else {
		cout << "error: smth went wrong((( \n file can not be open" << endl;
		return NULL;
	}
}

void print_Sys_Matrix(double**& A, size_t& n) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n + 1; ++j) {
			if (j == n)
				cout << setw(10) << "|";
			cout << setw(15) << A[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
}

void print_Matrix(double**& A, size_t& n) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) 
			cout << setw(15) << A[i][j];		
		cout << "\n";
	}
	cout << "\n";
}

void print_Vector(double*& vec, size_t& n)
{
	cout << "( ";
	for (size_t i = 0; i < n; ++i) 
		i != n-1 ? cout << vec[i] << ",\t" : cout << vec[i];
		
	cout << " )^T\n";
}



int copy_Matrix_A_to_B(double**& A, double**& B, size_t& n) {

	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
			B[i][j] = A[i][j];
	return 0;
}

int copy_Sys_A_to_B(double**& A, double**& B, size_t& n) {

	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n + 1; ++j)
			B[i][j] = A[i][j];
	return 0;
}

int copy_Vector_A_to_B(double*& A, double*& B, size_t& n) {

	for (size_t i = 0; i < n; ++i)		
			B[i] = A[i];
	return 0;
}


int marix_make_E(double** &A, size_t& n) {

	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
			(i == j) ? A[i][j] = 1.0 : A[i][j] = 0.0;

	return 0;
}

int matrix_mult(double** &A, double** &B, double**& T, size_t n)
{
	double temp = 0.0;

	//Transpose
	for (size_t i = 0; i < n; ++i)
		for (size_t j = i + 1; j < n; ++j)
		{
			temp = B[i][j];
			B[i][j] = B[j][i];
			B[j][i] = temp;
		}

	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
		{
			T[i][j] = 0.0;

			for (size_t k = 0; k < n; ++k)
				T[i][j] += A[i][k] * B[j][k];

			if (fabs(T[i][j]) < eps) {
				T[i][j] *= (1 / T[i][j] > 0) ? (1) : (-1);
			}
		}
	//������� ��������� ��� ������� B 
	//�������, ����������� �������� �� ����� ����� ����������, 
	//�� ��������� ������� ����������������, �� ��������� ������ ���
	for (size_t i = 0; i < n; ++i)
		for (size_t j = i + 1; j < n; ++j)
		{
			temp = B[i][j];
			B[i][j] = B[j][i];
			B[j][i] = temp;
		}


	return 0;
}

int matrix_vector_mult(double** A, double* B, double*& T, size_t n){


	for (size_t i = 0; i < n; ++i) {
		T[i] = 0.0;
		for (size_t j = 0; j < n; ++j) {

			T[i] += A[i][j] * B[j];
			if (fabs(T[i]) < eps) {
				T[i] *= (1 / T[i] > 0) ? (1) : (-1);
			}

		}
	}

	return 0;

}

int normalize_Sys_MatrixLine(size_t& k, double**& A, size_t& n) {

	double tmp = A[k][k];

	for (size_t j = 0; j < n + 1; ++j) {
		A[k][j] /= tmp;

		//solving -0 problem
		if (fabs(A[k][j]) < eps) {
			A[k][j] *= (1 / A[k][j] > 0) ? (1) : (-1);
		}
	}
	return 0;
}

int check_matrix_zero(double**& A, size_t& n) {
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
			if (fabs(A[i][j]) < eps) {
				A[i][j] = 0.0;
				A[i][j] *= (1 / A[i][j] > 0) ? (1) : (-1);
			}
	return 0;
}

int check_sys_matrix_zero(double**& A, size_t& n) {
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n + 1; ++j)
			if (fabs(A[i][j]) < eps && A[i][j] != 0.0) {
				A[i][j] = 0.0;
				A[i][j] *= (1 / A[i][j] > 0) ? (1) : (-1);
			}
	return 0;
}

int check_vector_zero(double*& A, size_t& n){

	for (size_t j = 0; j < n + 1; ++j)
		if (fabs(A[j]) < eps && A[j] != 0.0) {
			A[j] = 0.0;
			A[j] *= (1 / A[j] > 0) ? (1) : (-1);
		}

	return 0;
}

//--------//

int norm_1_vector(double*& v, double& norm, size_t& n) {

	norm = 0.0;
	for (size_t j = 0; j < n; ++j) {
		norm += fabs(v[j]);
	}

	if (fabs(norm) < eps && norm != 0.0) {
		norm = 0.0;
		norm *= (1 / norm > 0) ? (1) : (-1);
	}
	
	return 0;
}

int norm_inf_vector(double*& v, double& norm, size_t& n) {

	norm = 0.0;
	for (size_t j = 0; j < n; ++j) {
		if (fabs(fabs(norm) - fabs(v[j])) > eps) {
			if (fabs(norm) < fabs(v[j])) {
				norm = fabs(v[j]);
			}
		}
	}
	if (fabs(norm) < eps && norm != 0.0) {
		norm = 0.0;
		norm *= (1 / norm > 0) ? (1) : (-1);
	}
	return 0;
}

int norm_1_matrix(double**& A, double& norm, size_t& n) {

	norm = 0.0;
	double temp = 0.0;

	for (size_t i = 0; i < n; ++i) {
		temp = 0.0;
		for (size_t j = 0; j < n; ++j) {
			temp += fabs(A[j][i]);
		}		
		norm = (norm < temp) ? temp : norm;
	}
	if (fabs(norm) < eps && norm != 0.0) {
		norm = 0.0;
		norm *= (1 / norm > 0) ? (1) : (-1);
	}
	
	return 0;
}

int norm_inf_matrix(double**& A, double& norm, size_t& n) {

	norm = 0.0;
	double temp = 0.0;

	for (size_t i = 0; i < n; ++i) {
		temp = 0.0;
		for (size_t j = 0; j < n; ++j) {
			temp += fabs(A[i][j]);
		}
		norm = (norm < temp) ? temp : norm;
	}
	if (fabs(norm) < eps && norm != 0.0) {
		norm = 0.0;
		norm *= (1 / norm > 0) ? (1) : (-1);
	}

	return 0;
}

int discrepancy(double**& A, double*& x, double& d, size_t& n, int flag) {

	double* v = new double[n];
	matrix_vector_mult(A, x, v, n);
	//print_Sys_Matrix(A, n);
	for (size_t j = 0; j < n; ++j) {
		v[j] -= A[j][n];
		//cout << A[j][n] << "\n";
	}
	if (flag == 1) {
		norm_1_vector(v, d, n);
	}
	else
	{
		norm_inf_vector(v, d, n);
	}
	//print_Vector(x, n);
	//print_Vector(v, n);
	delete[] v;
	return 0;
}

int transpose_matrix(double**& A, size_t& n) {
	double temp;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			temp = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = temp;
		}
	}
	return 0;
}

int reverse_matrix(double**& A, double**& X, size_t& n) {

	double** B = new double* [n];
	for (size_t j = 0; j < n; ++j)
		B[j] = new double[n + 1];

	

	for (size_t j = 0; j < n; ++j)
	{
		copy_Matrix_A_to_B(A, B, n);

		for (size_t t = 0; t < n; ++t) {
			B[t][n] = 0.0;
		}
		B[j][n] = 1.0;
		
		method_gauss(B, n, X[j]);
	}
	
	transpose_matrix(X, n);
	check_matrix_zero(X, n);
	
	for (size_t j = 0; j < n; ++j)
		delete[] B[j];
	delete[] B;


	return 0;
}

int cond_matrix(double**& A, double& cond1, size_t& n, int flag) {

	double normA, normB;

	double** B = new double* [n];
	for (size_t j = 0; j < n; ++j)
		B[j] = new double[n+1];
	
	if (flag == 1) {
		norm_1_matrix(A, normA, n);
		
	}
	else
	{
		norm_inf_matrix(A, normA, n);
		
	}	
	
	reverse_matrix(A, B, n);
	
	if (flag == 1) {
		norm_1_matrix(B, normB, n);
	}
	else
	{
		norm_inf_matrix(B, normB, n);
	}
	
	cond1 = normA * normB;

	
	for (size_t j = 0; j < n; ++j)
		delete[] B[j];
	delete[] B;
	

	return 0;
}

int vector_A_min_B(double*& A, double* B, double*& res, size_t& n) {

	for (size_t j = 0; j < n; ++j) {
		res[j] = A[j] - B[j];
	}
	return 0;
}

int vector_Sys_A_min_B(double**& A, double** B, double*& res, size_t& n) {

	for (size_t j = 0; j < n; ++j) {
		res[j] = A[j][n] - B[j][n];
	}
	return 0;
}

int vector_valuation(double*& X, double*& _X, double& res, size_t& n, int flag) {

	double* vec = new double[n];
	double tmp = 0.0;

	vector_A_min_B(X, _X, vec, n);
	
	if (flag == 1) {
		norm_1_vector(vec, res, n);
		norm_1_vector(X, tmp, n);			
	}
	else {
		norm_inf_vector(vec, res, n);
		norm_inf_vector(X, tmp, n);
	}
	res /= tmp;

	delete[] vec;

	return 0;
}

int vector_valuation_Sys(double**& X, double**& _X, double& res, size_t& n, int flag) {

	double* vec = new double[n];
	double tmp = 0.0;


	vector_Sys_A_min_B(X, _X, vec, n);

	if (flag == 1) {
		norm_1_vector(vec, res, n);

		for (size_t i = 0; i < n; ++i)
			vec[i] = X[i][n];

		norm_1_vector(vec, tmp, n);
		
	}
	else {
		norm_inf_vector(vec, res, n);

		for (size_t i = 0; i < n; ++i)
			vec[i] = X[i][n];

		norm_inf_vector(vec, tmp, n);
		
	}
	res /= tmp;

	delete[] vec;

	return 0;
}

int vozm(double**& A, double*& x_save, double & max1, double & max_inf, size_t& n) {
	

	double pl = 0.01;
	max1 = 0.0;
	max_inf = 0.0;

	double temp1 = 0.0, temp_inf = 0.0;


	double** A_save = new double* [n];
	for (size_t j = 0; j < n; ++j)
		A_save[j] = new double[n + 1];

	double delta_x = 0.0, delta_b = 0.0;

	double* x = new double [n];

	cout << "1\n";
	for (size_t l = 0; l < 10; ++l) {

		srand(time(NULL));

		copy_Sys_A_to_B(A, A_save, n);

		double delta_b =0.0;
		//cout << "2\n";
		for (size_t i = 0; i < n; ++i) {
			//cout << i <<"aafdfdfaa\n";
			delta_b = pl * pow(10 + l, rand() % -1 + 1);
			A_save[i][n] += delta_b;

		}
		//cout << "3-1\n";
		method_gauss(A_save, n, x);
		//cout << "3-2\n";

		vector_valuation(x_save, x, delta_x, n, 1);
		vector_valuation_Sys(A_save, A, delta_b, n, 1);
		//cout << "4\n";
		temp1 = delta_x / delta_b;

		//cout << temp1 << "  1  " << max1 << "\n";

		if (fabs(temp1 - max1) > eps) {		
			
			if (temp1 > max1) {
				max1 = temp1;
			}
		}
		//cout << "5\n";
		vector_valuation(x_save, x, delta_x, n, 0);
		vector_valuation_Sys(A_save, A, delta_b, n, 0);
		//cout << "6\n";

		temp_inf = delta_x / delta_b;

		//cout << temp_inf << "  inf  " << max_inf << "\n";
		//cout << "7\n";
		if (fabs(temp_inf - max_inf) > eps) {
			
			if (temp_inf > max_inf) {
				max_inf = temp_inf;
			}

		}
		

	}

	return 0;

}