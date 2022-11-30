#include <iostream>
#include <clocale>
#include <cmath>

#include "Matrix_public.h"
#include "interpol.h"
#include "spline.h"

#define pi 3.14159
using namespace std;

double func(double& x){
	return exp(x) + x*x + 35;
}

int main(int args, char** argv){

	setlocale(LC_ALL, "Rus");

	size_t n = 10; //vector size
	
	
	double* x = new double[n];
	double xx = 0.222;

	cout << "x = " << xx << "\n";

	
	double a = -1.0, b = 1.0;

	//ravn

	double h = (b - a) / (n-1);
	for(size_t i = 0; i <= n; ++i){
		x[i] = a + i * h;
	}

	double* y = new double[n];
	for(size_t i = 0; i < n; ++i){
		y[i] = func(x[i]);
	}
	cout << "Lagrange interpolation (unoform grid): \n";
	double k = interpol_lagrange(xx, x ,y, n);
	cout << k << endl;


	//cheb

	for(size_t i = 0; i <= n; ++i){
		x[i] = (a + b) / n + (b - a) / n * cos((2*i + 1) * pi / 2 / (n+1));
	}

	for(size_t i = 0; i < n; ++i){
		y[i] = func(x[i]);
	}

	cout << "Lagrange interpolation (Chebyshov's grid): \n";
	k = interpol_lagrange(xx, x ,y, n);
	cout << k << endl;


	string path = "/home/san/Code/labs_comput/lab_3/test_files/interpol.txt";

	double* x_new = read_vec(path, n);
	for(size_t i = 0; i < n; ++i){
		y[i] = func(x_new[i]);
	}

	cout << "Extrapolation: \n";

	k = interpol_lagrange(xx, x_new ,y, n);
	cout << k << endl;


	h = (b - a) / (n-1);
	for(size_t i = 0; i < n; ++i){
		x[i] = a + i * h;
	}
	
	for(size_t i = 0; i < n; ++i){
		y[i] = func(x[i]);
	}

	double** A = spline(xx, x, y, n);

	int kk = -1;
    for(size_t i = 0; i < n-1; ++i){
        if(x[i+1] - xx > 0 && xx - x[i] > 0){ 
			double res = A[i][0] + A[i][1] * (xx - x[i]) + A[i][2] * pow(xx - x[i], 2) + A[i][3] * pow(xx - x[i], 3);
			cout << "Spline (uniform grid): \n";			
			cout << res << "\n";
            break;
		}
		
    }


	// for(size_t i = 0; i < n-1; ++i){
	// 	x[i] = (a + b) / n + (b - a) / n * cos((2*i + 1) * pi / 2 / (n+1));
	// }

	// print_Vector(x, n);

	// for(size_t i = 0; i < n; ++i){
	// 	y[i] = func(x[i]);
	// }
	
	// A = spline(xx, x, y, n);

	// kk = -1;
    // for(size_t i = 0; i < n-1; ++i){
    //     if(x[i+1] - xx > 0){ 
	// 		double res = A[i][0] + A[i][1] * (xx - x[i]) + A[i][2] * pow(xx - x[i], 2) + A[i][3] * pow(xx - x[i], 3);
	// 		cout << "Spline (Chebyshov's grid): \n";			
	// 		cout << res << "\n";
    //         break;
	// 	}
		
    // }

		
	cout << "--------------------------------\n";
	cout << "--------------------------------\n";

	return 0;
}