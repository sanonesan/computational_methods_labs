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

	double** A = readMatrixFromFile(path, n); 
		
	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "System:\n";
	print_Matrix(A, n);

	double* lambda_mas = hessenberg(A, n);

	for(size_t i = 0; i < n; ++i){
        cout << "Lambda_" << i << " = " << lambda_mas[i] << "\n";
    }

	double* x0 = new double [n];

	for(size_t i = 0; i < n; ++i){
		x0[i] = 0;
	}
	x0[2]=1;


	double* x; 
	cout << "-----------------------------------------------------\n";
	cout << "\ninverse iteration: \n\n";
	for(size_t i = 0; i < n; ++i){
		cout << "Lambda_" << i << " = " << lambda_mas[i] << "\n";
		x = inverse_iter_method(A, x0, lambda_mas[i], n);
		print_Vector(x, n);
		cout << "\n";
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
	}

	delete[] A;
	delete[] x0;
	delete[] x;
	
	return 0;
}
