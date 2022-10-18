#include <iostream>

#include "Matrix_public.h"
#include "method_gauss.h"
#include "simple_iter_method.h"
#include "Jacoby_method.h"
#include <clocale>

using namespace std;



int main(int args, char** argv){

	setlocale(LC_ALL, "Rus");


	size_t n; //vector size
	//string path = "C:\\Code\\labs_comput\\lab_1\\input.txt";
	//string path_dis = "C:\\Code\\labs_comput\\lab_1\\input_disturbed.txt";

	string path = "/home/san/Code/labs_comput/lab_2/test_files/data35.txt";
	string path_dis = "/home/san/Code/labs_comput/lab_2/test_files/P_DATA6_D.TXT";

	double** A = readMatrixFromFile(path, n); //Matrix A|b of sys Ax=b
	
	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "System:\n";
	double* x0 = new double [n];
	// print_Sys_Matrix(A, n);

	// method_gauss(A, n, x0);
	// print_Vector(x0, n);
	
	print_Sys_Matrix(A, n);

	A = readMatrixFromFile(path, n); //Matrix A|b of sys Ax=b
	
	double norm = 0.0;
	norm_1_matrix(A, norm, n);
	cout << norm << "\n" ;
	double* x;
	//double* x0 = new double [n];
	
	for (size_t i = 0; i < n; ++i){
		x0[i] = 1;
	}
	

	x = simple_iter_method(A, x0, 0.051282, n);

	print_Vector(x, n);

	x = Jacoby_method(A, x0, n);

	print_Vector(x, n);
	
	// double** L;
	// double** D;
	// double** U;

	//LDU_method(A, L, D, U, n);


	// reverse_matrix(D, L, n);
	// print_Matrix(L, n);
	// reverse_matrix_qr(D, L, n);
	// print_Matrix(L, n);



	for (size_t j = 0; j < n; ++j) {
		delete[] A[j];
	}

	delete[] A;
	
	return 0;
}
