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
	//string path = "C:\\Code\\labs_comput\\lab_1\\input.txt";
	//string path_dis = "C:\\Code\\labs_comput\\lab_1\\input_disturbed.txt";

	string path = "/home/san/Code/labs_comput/lab_3/test_files/lab_4_mat.txt";
	//string path = "/home/san/Code/labs_comput/lab_3/test_files/test311.txt";

	//string path_dis = "/home/san/Code/labs_comput/lab_2/test_files/P_DATA6_D.TXT";

	double** A = readMatrixFromFile(path, n); //Matrix A|b of sys Ax=b
	
	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	cout << "System:\n";
	print_Matrix(A, n);

	//double g = 0.5;
	double* lambda_mas = hessenberg(A, n);

	for(size_t i = 0; i < n; ++i){
        cout << "Lambda_" << i << " = " << lambda_mas[i] << "\n";
    }

	double* x0 = new double [n];

	for(size_t i = 0; i < n; ++i){
		x0[i] = 1;
	}


	double* x = inverse_iter_method(A, x0, lambda_mas[0], n);

	print_Vector(x, n);


	// x0[0] = 0.8;
	// x0[1] = 0.0;
	// x0[2] = 0.2;
	// x0[3] = 0.4;


	// double* x ;
	
	// for(size_t i = 0; i < 1; ++i){
	// 	x = inverse_iter_method_Rayleigh(A, x0, lambda_mas[0], n);
	// 	print_Vector(x, n);
	// }


	for (size_t j = 0; j < n; ++j) {
		delete[] A[j];
	}

	delete[] A;
	delete[] x0;
	delete[] x;
	
	return 0;
}
