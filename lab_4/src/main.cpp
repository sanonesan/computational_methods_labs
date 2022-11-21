#include <iostream>

#include "Matrix_public.h"
#include "eigenValues_QR_hessenberg.h"
#include "inverse_iter_method.h"
#include <clocale>
//#include <time.h>


using namespace std;

// double inaccuracy(double**& A, double*& x){
// 	double q = 0.0;

// 	for 
// }

int main(int args, char** argv){

	setlocale(LC_ALL, "Rus");


	size_t n; //vector size

	string path = "/home/san/Code/labs_comput/lab_4/test_files/lab_4_mat.txt";
	//string path = "/home/san/Code/labs_comput/lab_4/test_files/input1.txt";


	double** A = readMatrixFromFile(path, n); 
		
	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "System:\n";
	print_Matrix(A, n);


	double** A_save = new double* [n];
	for(size_t i = 0; i < n; ++i){
		A_save[i] = new double [n];
	}

	//copy_Matrix_A_to_B(A, A_save, n);

	//transpose_matrix(A, n);

	double* lambda_mas = hessenberg(A, n);

	for(size_t i = 0; i < n; ++i){
        cout << "Lambda_" << i << " = " << lambda_mas[i] << "\n";
    }

	double* x0 = new double [n];

	for(size_t i = 0; i < n; ++i){
		x0[i] = 0;
	}
	
	x0[0]=0.8;
	x0[1]=0;
	x0[2]=0.2;
	x0[3]=0.4;


// Lambda_0 = 0.997313
// iter = 3
// oper = 122
// ( 0.864461,	0.00338916,	0.244595,	0.439169 )^T

// Lambda_1 = 2.00425
// iter = 5
// oper = 200
// ( 0.0119071,	-0.709759,	0.607555,	-0.356337 )^T

// Lambda_2 = 2.98702
// iter = 5
// oper = 200
// ( -0.502548,	-0.0157964,	0.430896,	0.749349 )^T

// Lambda_3 = 4.01142
// iter = 2
// oper = 80
// ( -0.00343211,	0.704259,	0.620789,	-0.344426 )^T

	double* x; 
	double* check = new double [n];
	cout << "-----------------------------------------------------\n";
	cout << "\ninverse iteration: \n\n";
	for(size_t i = 0; i < n; ++i){
		cout << "Lambda_" << i << " = " << lambda_mas[i] << "\n";
		x = inverse_iter_method(A, x0, lambda_mas[i], n);

		//cout << "\nEigenvector (A^T) * x' = lambda * x': \n";
		print_Vector(x, n);
		
		// cout << "\nx' * A: \n";
		// vector_matrix_mult(x, A_save, check, n);
		// print_Vector(check,n);

		// cout << "x' * lambda: \n";

		// for(size_t j = 0; j < n; ++j){
		// 	cout << x[j] * lambda_mas[i] << " ";
		// }
		cout << "\n\n";
	}
	
	cout << "-----------------------------------------------------\n";
	cout << "\nRayleigh: \n\n";
	for(size_t i = 0; i < n; ++i){
		cout << "Lambda_" << i << " = " << lambda_mas[i] << "\n";
		x = inverse_iter_method_Rayleigh(A, x0, lambda_mas[i], n);
		print_Vector(x, n);
		cout << "\n";
	}


	for (size_t j = 0; j < n; ++j) {
		delete[] A[j];
		delete[] A_save[j];
	}

	delete[] A;
	delete[] A_save;
	delete[] check;
	delete[] x0;
	delete[] x;
	
	return 0;
}
