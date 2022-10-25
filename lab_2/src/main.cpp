#include <iostream>

#include "Matrix_public.h"
#include "method_gauss.h"
#include "simple_iter_method.h"
#include "Jacoby_method.h"
#include "Relaxation_method.h"
#include <clocale>
#include <time.h>


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
	
	//double norm = 0.0;
	//norm_1_matrix(A, norm, n);
	//cout << norm << "\n" ;
	double* x;
	//double* x0 = new double [n];
	
	for (size_t i = 0; i < n; ++i){
		x0[i] = 1;
	}
	
	cout << "\n-----------Simple iter------------\n";
	
	cout << "tau = 0.051282\n";
	x = simple_iter_method(A, x0, 0.051282, n);

	print_Vector(x, n);

	cout << "\n-----------Simple iter------------\n";
	
	cout << "tau = 0.052\n";
	x = simple_iter_method(A, x0, 0.052, n);

	print_Vector(x, n);

	cout << "\n--------------Jacoby---------------\n";


	x = Jacoby_method(A, x0, n);

	print_Vector(x, n);



	cout << "\n--------------Relaxation---------------\n";

	double w = 0.7;

	cout << "omega = " << w << "\n";
	

	x = Relaxation_method(A, x0, w, n);

	print_Vector(x, n);

	cout << "\n--------------Relaxation---------------\n";

	w = 0.2;

	cout << "omega = " << w << "\n";
	

	x = Relaxation_method(A, x0, w, n);

	print_Vector(x, n);




	cout << "\n--------------Seidel---------------\n";

	w = 1;

	//cout << "omega = " << w << "\n";

	x = Relaxation_method(A, x0, w, n);

	print_Vector(x, n);


	cout << "\n----------3 diag ------------\n";

	size_t k = 250;
	double* a = new double[k];
	double* b = new double[k];
	double* c = new double[k];
	double* d = new double[k];

	for(size_t i = 0; i < k; ++i){
		a[i] = 6;
		b[i] = 8;
		c[i] = 1;
		d[i] = i+1;
	}

	clock_t start = clock();

	x = Relaxation_method(a, b, c, d, x0, w, k);

	clock_t end = clock();

	double seconds = (double)(end - start) / CLOCKS_PER_SEC;

	//cout << "The time: " << seconds << " \n ------------------------- \n";
	//print_Vector(x, k);
	
	cout << "\n----------3 diag transp------------\n";

	for(size_t i = 0; i < k; ++i){
		a[i] = 1;
		b[i] = 8;
		c[i] = 6;
		d[i] = i+1;
	}
	start = clock();
	x = Relaxation_method(a, b, c, d, x0, w, k);
	end = clock();

	seconds = (double)(end - start) / CLOCKS_PER_SEC;
	//cout << "The time: " << seconds << " \n";
	cout << "\n-------------------------------\n";
	//print_Vector(x, k);

	for (size_t j = 0; j < n; ++j) {
		delete[] A[j];
	}

	delete[] A;
	
	return 0;
}
