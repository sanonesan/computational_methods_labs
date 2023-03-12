#include <iostream>

#include "Matrix_public.h"
#include "method_gauss.h"
#include "QR_method.h"
#include <clocale>
#include <cmath>

using namespace std;

int main(int args, char** argv){

	setlocale(LC_ALL, "Rus");

	size_t n; 

	string path = "../test_files/P_DATA6.TXT";
	string path_dis = "../test_files/P_DATA6_D.TXT";

	double** A = readMatrixFromFile(path, n); //Matrix A|b of sys Ax=b
	
	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "System:\n";

	print_Sys_Matrix(A, n);

	double* x = new double [n]; //Answer for: Ax=b or QRx=b
	double* x_save = new double[n];

	cout << "\nGauss-method:\n";


	int check = method_gauss(A, n, x);


	if (check == 0) {
		print_Sys_Matrix(A, n);
		cout << "Result: \n";
		print_Vector(x, n);		

		for (size_t i = 0; i < n; ++i)
			x_save[i] = x[i];
	}


	double** Q = new double* [n];
	double** R = new double* [n];
	for (size_t j = 0; j < n; ++j) {
		Q[j] = new double[n];
		R[j] = new double[n];
	}
		

	A = readMatrixFromFile(path, n);
	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "\nSystem:\n";

	print_Sys_Matrix(A,n);

	cout << "\nQR-method:\n";
	check = QR_method_find_Q_and_R(
		Q, 
		R, 
		A, n, x);
	if (check == 0) {
		print_Sys_Matrix(A, n);
		cout << "Q: \n";
		print_Matrix(Q, n);
		cout << "R: \n";
		print_Matrix(R, n);

		cout << "Result: \n";
		print_Vector(x, n);
	}

	A = readMatrixFromFile(path, n);

	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "\n\nЧисленное решение: \n";
	print_Vector(x, n);
	cout << "\n Невязка: \n";

	double dis_1, dis_inf; 
	discrepancy(A, x, dis_1, n, 1);
	discrepancy(A, x, dis_inf, n, 0);

	cout << "  norm_1: " << dis_1 << "\n" << "norm_inf: " << dis_inf << "\n\n";
	cout << "Точное значение числа обусловленности: \n";

	double cond1 = 0.0, condinf = 0.0;
	cond_matrix(A, cond1, n, 1);
	//cout << "a: \n";
	cond_matrix(A, condinf, n, 0);

	double max1 = 0.0, max_inf = 0.0;
	vozm(A, x, max1, max_inf, n);

	//cout << max1
	//vector_valuation(x_save, x, delta_x, n, 1);
	//vector_valuation_Sys(A_save, A, delta_b, n, 1);

	cout << "vozm: \n";
	cout << "  norm_1: delta_x / delta_b = " << max1 << endl;
	cout << "norm_inf: delta_x / delta_b = " << max_inf << endl;

	//cout << "a: \n";

	cout << "  norm_1: " << cond1 << "\n" << "norm_inf: " << condinf << "\n";
	
	cout << "\n\nРешение системы с возмущенной правой частью: \n";

	

	A = readMatrixFromFile(path_dis, n); //Matrix A|b of sys Ax=b

	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "System:\n";

	print_Sys_Matrix(A, n);

	


	cout << "\nGauss-method:\n";


	check = method_gauss(A, n, x);


	if (check == 0) {
		print_Sys_Matrix(A, n);
		cout << "Result: \n";
		print_Vector(x, n);
	}


	double** A_save = readMatrixFromFile(path, n);

	A = readMatrixFromFile(path_dis, n);
	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "\nSystem:\n";

	print_Sys_Matrix(A, n);

	cout << "\nQR-method:\n";
	check = QR_method_find_Q_and_R(
		Q,
		R,
		A, n, x);
	if (check == 0) {
		print_Sys_Matrix(A, n);
		cout << "Q: \n";
		print_Matrix(Q, n);
		cout << "R: \n";
		print_Matrix(R, n);

		cout << "Result: \n";
		print_Vector(x, n);
	}

	A = readMatrixFromFile(path_dis, n);

	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "--------------------------------\n";

	cout << "\n\nЧисленное решение: \n";
	print_Vector(x, n);
	cout << "\n Невязка: \n";
	
	discrepancy(A, x, dis_1, n, 1);
	discrepancy(A, x, dis_inf, n, 0);
	cout << "  norm_1: " << dis_1 << "\n" << "norm_inf: " << dis_inf << "\n\n";
	cout << "Значение числа обусловленности: \n";
	
	cond_matrix(A, cond1, n, 1);
	cond_matrix(A, condinf, n, 0);
	
	cout << "  norm_1: " << cond1 << "\n" << "norm_inf: " << condinf << "\n\n";


	double delta_x = 0.0, delta_b = 0.0;

	vector_valuation(x_save, x, delta_x, n, 1);
	vector_valuation_Sys(A_save, A, delta_b, n, 1);
	vozm(A, x, max1, max_inf, n);
	
	cout << "Значение числа обусловленности: \n";
	cout << "vozm: \n";

	if (max1 < delta_x / delta_b)
		cout << "  norm_1: delta_x / delta_b = " << delta_x << " / " << delta_b  << " = " << delta_x / delta_b << endl;
	else
		cout << "  norm_1: delta_x / delta_b = " << max1 << "\n";


	vector_valuation(x_save, x, delta_x, n, 0);
	vector_valuation_Sys(A_save, A, delta_b, n, 0);


	if (max_inf < delta_x / delta_b)
		cout << "norm_inf: delta_x / delta_b = " << delta_x << " / " << delta_b << " = " << delta_x / delta_b << endl;
	else
		cout << "norm_inf: delta_x / delta_b = " << max_inf << endl;


	//double max1 = 0.0, max_inf = 0.0;

	double** B = new double* [n];
	for (size_t j = 0; j < n; ++j)
		B[j] = new double[n];

	cout << "\n-----------------------\n";
	reverse_matrix(A, B, n);
	print_Matrix(B, n);
	cout << "\n-----------------------\n";
	reverse_matrix_qr(A, B, n);
	print_Matrix(B, n);


	for (size_t j = 0; j < n; ++j)
		delete[] B[j];
	delete[] B;

	
	//cout << max1
	//vector_valuation(x_save, x, delta_x, n, 1);
	//vector_valuation_Sys(A_save, A, delta_b, n, 1);


	//vector_valuation(x_save, x, delta_x, n, 0);
	//vector_valuation_Sys(A_save, A, delta_b, n, 0);


	for (size_t j = 0; j < n; ++j) {
		delete[] A[j];
		delete[] A_save[j];
	}

	delete[] A;
	delete[] A_save;
	delete[] x;
	delete[] x_save;



	return 0;
}
