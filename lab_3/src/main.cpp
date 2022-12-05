#include <iostream>
#include <clocale>
#include <cmath>

#include "Matrix_public.h"
#include "interpol.h"
#include "spline.h"

#include <fstream>
#include <vector>


#define pi 3.14159
using namespace std;

double func(double& x){
	//return exp(x) + x*x + 35;
	return exp(x);
}

double func1(double& x){
	//return exp(x) + x*x + 35;
	return x;
}

double func2(double& x){
	//return exp(x) + x*x + 35;
	return x * x;
}

double func3(double& x){
	//return exp(x) + x*x + 35;
	return 1 / (1 + x * x);
}

int main(int args, char** argv){

	setlocale(LC_ALL, "Rus");
	ofstream file;
	string path1;

	size_t n = 10; //vector size
	
	
	double xx = 2.2;

	cout << "x = " << xx << "\n";

	
	// double a = 0.0, b = 2.0;
	double a = -5.0, b = 5.0;


	//ravn

	size_t c_size = 1000;
	double* cords_interpol = new double [c_size];
	double h = ((b) - (a)) / (c_size-1);
	for(size_t i = 0; i < c_size; ++i){
		cords_interpol[i] = (a) + i * h;
	}

	h = (b - a) / (n-1);
	// h = 0.2;
	// n = (b-a) / h + 1;
	double* x = new double[n];

	for(size_t i = 0; i < n; ++i){
		x[i] = a + i * h;
	}

	double* y = new double[n];
	for(size_t i = 0; i < n; ++i){
		y[i] = func3(x[i]);
	}

/*-----------------------

	UNIFORM GRID

------------------------*/

	

	

	// //print_Vector(x, n);

	
	cout << "Lagrange interpolation (uniform grid): \n";
	double k = interpol_lagrange(xx, x ,y, n);
	cout << "L(2.2) = " << k << "\n";
	cout << "exp(2.2) = " << func(xx) << "\n";
	cout << "| exp(2.2) - L(2.2) | = " << fabs(k - func(xx))  << "\n";




	path1 = "/home/san/Code/labs_comput/lab_3/interpol_uniform_grid_out.txt";

    file.open(path1);


    if (file.is_open())
    {	
		
		for(size_t i = 0; i < c_size; ++i){			
			file << cords_interpol[i] << " " << interpol_lagrange(cords_interpol[i], x ,y, n) << "\n";
		}
		file.close();
    }

/*-----------------------

	CHEBYSHOV GRID

------------------------*/

	for(size_t i = 0; i <= n; ++i){
		x[i] = (a + b) / n + (b - a) / n * cos((2*i + 1) * pi / 2 / (n+1));
	}

	for(size_t i = 0; i < n; ++i){
		y[i] = func3(x[i]);
	}

	cout << "Lagrange interpolation (Chebyshov's grid): \n";
	k = interpol_lagrange(xx, x ,y, n);
	cout << k << endl;

	path1 = "/home/san/Code/labs_comput/lab_3/interpol_chebyshov_grid_out.txt";

    file.open(path1);

    if (file.is_open())
    {	
		
		for(size_t i = 0; i < c_size; ++i){			
			file << cords_interpol[i] << " " << interpol_lagrange(cords_interpol[i], x ,y, n) << "\n";
		}
		file.close();
    }



/*-----------------------

	EXTRAPOLATION

------------------------*/

	string path = "/home/san/Code/labs_comput/lab_3/test_files/interpol.txt";

	double* x_new = read_vec(path, n);
	for(size_t i = 0; i < n; ++i){
		y[i] = func3(x_new[i]);
	}

	cout << "Extrapolation: \n";

	k = interpol_lagrange(xx, x_new ,y, n);
	cout << k << endl;


	h = (b - a) / (n-1);
	//h = 0.2;

	for(size_t i = 0; i < n; ++i){
		x[i] = a + i * h;
	}
	
	for(size_t i = 0; i < n; ++i){
		y[i] = func3(x[i]);
	}

	path1 = "/home/san/Code/labs_comput/lab_3/extrapol_out.txt";

    file.open(path1);

    if (file.is_open())
    {	
		
		for(size_t i = 0; i < c_size; ++i){			
			file << cords_interpol[i] << " " << interpol_lagrange(cords_interpol[i], x ,y, n) << "\n";
		}
		file.close();
    }


/* SPLINE */

	double** A = spline(x, y, n);
	xx = 0.4597;
	int kk = -1;
	double res = 0.0;
    for(size_t i = 0; i < n-1; ++i){
        if(x[i+1] - xx > 0 && xx - x[i] > 0){ 
			res = A[i][0] + A[i][1] * (xx - x[i]) + A[i][2] * pow(xx - x[i], 2) + A[i][3] * pow(xx - x[i], 3);
			cout << "Spline (uniform grid): \n";			
			cout << fabs(res - xx) << "\n";
			cout << fabs(res - xx*xx) << "\n";

            break;
		}
		
    }

	path1 = "/home/san/Code/labs_comput/lab_3/spline_out.txt";

    file.open(path1);

    if (file.is_open())
    {	
		h = ((b) - (a)) / (c_size-1);
		for(size_t i = 0; i < c_size; ++i){
			cords_interpol[i] = (a) + i * h;
		}

		
		
		for(size_t j = 0; j < c_size; ++j){	

			xx = cords_interpol[j];
			for(size_t i = 0; i < n-1; ++i){
				
				if(x[i+1] - xx > 0 && xx - x[i] > 0){ 
					file << xx << " " << A[i][0] + A[i][1] * (xx - x[i]) + A[i][2] * pow(xx - x[i], 2) + A[i][3] * pow(xx - x[i], 3) << "\n";
					break;
				}
		
			}
		}
		file.close();
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
