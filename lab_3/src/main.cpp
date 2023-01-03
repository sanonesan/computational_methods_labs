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

double func(double &x)
{
	// return exp(x) + x*x + 35;
	//return exp(x);
	return 1 / (atan(1 +10* x*x));
}

double func1(double &x)
{
	// return exp(x) + x*x + 35;
	return x;
}

double func2(double &x)
{
	// return exp(x) + x*x + 35;
	return x * x;
}

double func3(double &x)
{
	// return exp(x) + x*x + 35;
	return 1 / (1 + x * x);
}

double func_pow2(double &x, int k)
{
	return pow(pow(2, k), x);
}

int main(int args, char **argv)
{

	setlocale(LC_ALL, "Rus");
	ofstream file;
	string path1, path2;

	size_t n = 10; // vector size
	double xx = 1;
	// double a = 0.0, b = 2.0;

	double a = -5.0, b = 5.0;

	size_t c_size = 1000;
	double *cords_interpol = new double[c_size];
	double h = ((b) - (a)) / (c_size - 1);
	for (size_t i = 0; i < c_size; ++i)
	{
		cords_interpol[i] = (a) + i * h;
	}

	h = (b - a) / (n - 1);
	// h = 0.2;
	// n = (b-a) / h + 1;
	double *x = new double[n];

	for (size_t i = 0; i < n; ++i)
	{
		x[i] = a + i * h;
	}

	double *y = new double[n];
	for (size_t i = 0; i < n; ++i)
	{
		y[i] = func3(x[i]);
	}

	path1 = "/home/san/Code/labs_comput/lab_3/interpol_uniform_grid_out.txt";
	path2 = "/home/san/Code/labs_comput/lab_3/interpol_chebyshov_grid_out.txt";

	double k = 0;

	xx = 0.348;
	cout << "x = " << xx << "\n";

	/*-----------------------

		UNIFORM GRID

	------------------------*/

	size_t u = 8;
	cout << "Lagrange interpolation (uniform grid): \n";

	double *mas = new double[u];

	for (size_t m = 1; m <= u; ++m)
	{

		n = pow(2, m)+1;

		delete[] x;
		delete[] y;

		x = new double[n];
		y = new double[n];

		h = (b - a) / (n-1);

		for (size_t i = 0; i < n; ++i)
		{
			x[i] = a + i * h;
		}

		for (size_t i = 0; i < n; ++i)
		{
			y[i] = func(x[i]);
		}
		mas[m - 1] = fabs(func(xx) - interpol_lagrange(xx, x, y, n));
	}
	cout << "Вектор нормы ошибки: \n";
	print_Vector(mas, u);

	///////////////////////////////////////////////////


	// cout << "L(2.2) = " << k << "\n";
	// cout << "exp(2.2) = " << func(xx) << "\n";
	// cout << "| exp(2.2) - L(2.2) | = " << fabs(k - func(xx))  << "\n";

	// file.open(path1);

	// if (file.is_open())
	// {
	// 	for(size_t i = 0; i < c_size; ++i){
	// 		file << cords_interpol[i] << " " << interpol_lagrange(cords_interpol[i], x ,y, n) << "\n";
	// 	}
	// 	file.close();
	// }

	/*-----------------------

		CHEBYSHOV GRID

	------------------------*/

	cout << "\nLagrange interpolation (Chebyshov's grid): \n";

	for (size_t m = 1; m <= u; ++m)
	{

		n = pow(2, m);
		// delete[] x;
		// delete[] y;

		x = new double[n];
		y = new double[n];

		for (size_t i = 0; i <= n; ++i)
		{
			x[i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * pi / 2 / (n + 1));
		}

		for (size_t i = 0; i < n; ++i)
		{
			y[i] = func(x[i]);
		}
		mas[m - 1] = fabs(func(xx) - interpol_lagrange(xx, x, y, n));
	}
	cout << "Вектор нормы ошибки: \n";
	print_Vector(mas, u);

	n = 10;

	x = new double[n];
	y = new double[n];

	for(size_t i = 0; i < n; ++i){
		x[i] = (a + b) / 2 + (b - a) / 2 * cos((2*i + 1) * pi / 2 / (n));
	}

	//print_Vector(x, n);

	for(size_t i = 0; i < n; ++i){
		//y[i] = 1;
		y[i] = func3(x[i]);
	}

	cout << "Lagrange interpolation (Chebyshov's grid): \n";
	k = interpol_lagrange(xx, x ,y, n);
	cout << k << endl;

	file.open(path2);

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

	// string path = "/home/san/Code/labs_comput/lab_3/test_files/interpol.txt";

	// double* x_new = read_vec(path, n);
	// for(size_t i = 0; i < n; ++i){
	// 	y[i] = func3(x_new[i]);
	// }

	// cout << "Extrapolation: \n";

	// k = interpol_lagrange(xx, x_new ,y, n);
	// cout << k << endl;

	// h = (b - a) / (n-1);
	// //h = 0.2;

	// for(size_t i = 0; i < n; ++i){
	// 	x[i] = a + i * h;
	// }

	// for(size_t i = 0; i < n; ++i){
	// 	y[i] = func3(x[i]);
	// }

	// path1 = "/home/san/Code/labs_comput/lab_3/extrapol_out.txt";

	// file.open(path1);

	// if (file.is_open())
	// {

	// 	for(size_t i = 0; i < c_size; ++i){
	// 		file << cords_interpol[i] << " " << interpol_lagrange(cords_interpol[i], x ,y, n) << "\n";
	// 	}
	// 	file.close();
	// }

	/* SPLINE */

	//Оценка P

	n = 10;

	h = (b - a) / (n - 1);

	delete[] x;
	delete[] y;

	x = new double[n];
	y = new double[n];

	for (size_t i = 0; i < n; ++i)
	{
		x[i] = a + i * h;
	}

	for (size_t i = 0; i < n; ++i)
	{
		y[i] = func2(x[i]);
	}

	double **A = spline(x, y, n);

	xx = 0.4597;

	int kk = -1;
	double res = 0.0;
	for (size_t i = 0; i < n - 1; ++i)
	{
		if (x[i + 1] - xx > 0 && xx - x[i] > 0)
		{
			res = A[i][0] + A[i][1] * (xx - x[i]) + A[i][2] * pow(xx - x[i], 2) + A[i][3] * pow(xx - x[i], 3);
			cout << "\nSpline (uniform grid): \n";
			cout << res << "\n adgdsagasdg \n";

			cout << fabs(res - xx*xx) << "\n";
			break;
		}
	}

	n = 16;
	delete[] x;
	delete[] y;
	x = new double[n];
	y = new double[n];

	for (size_t i = 0; i < n; ++i)
	{
		x[i] = a + i * h;
	}

	for (size_t i = 0; i < n; ++i)
	{
		y[i] = func1(x[i]);
	}

	A = spline(x, y, n);

	double res1 = 0.0;
	for (size_t i = 0; i < n - 1; ++i)
	{
		if (x[i + 1] - xx > 0 && xx - x[i] > 0)
		{
			res1 = A[i][0] + A[i][1] * (xx - x[i]) + A[i][2] * pow(xx - x[i], 2) + A[i][3] * pow(xx - x[i], 3);
			cout << "Spline (uniform grid): \n";
			cout << res1 << "\n";
			cout << res1 << "\n";
			cout << fabs(res - xx) << "\n";
			// cout << fabs(res - xx*xx) << "\n";
			break;
		}
	}

	cout << fabs(res1 - res) << endl;

	cout << "\n" << log2(fabs(func(xx) - res) / fabs(func(xx) - res1)) << "\n";

	path1 = "/home/san/Code/labs_comput/lab_3/spline_out.txt";

	file.open(path1);

	for (size_t i = 0; i < n; ++i)
	{
		y[i] = func2(x[i]);
	}

	A = spline(x, y, n);

	cout << "\n\n\n\nDZ2: \n";
	for(size_t i = 0; i < n-1; ++i){
		for(size_t j = 0; j < 4; ++j){
			cout << A[i][j] << "    ";
		}
		cout << "\n";
	}


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

	cout << "--------------------------------\n";
	cout << "--------------------------------\n";
	xx = 2;
	// cout << func_pow2(xx, 2) << "\n";

	n = 4;
	delete[] x;
	delete[] y;
	x = new double[n];
	y = new double[n];

	x[0] = 0;
	x[1] = 0.333333;
	x[2] = 0.666667;
	x[3] = 1;

	y[0] = 0.0;
	y[1] = -0.000843679;
	y[2] = -0.00325299;
	y[3] = -0.00712056;

	print_Vector(x, n);
	print_Vector(y,n);

	A = spline(x, y, n);

	cout << "\n\n\n\nDZ2: \n";
	for(size_t j = 0; j < 4; ++j){
		if (j == 0)
			cout << "a: ";
		if (j == 1)
			cout << "b: ";
		if (j == 2)
			cout << "c: ";
		if (j == 3)
			cout << "d: ";
		for(size_t i = 0; i < n-1; ++i){
			
			cout << A[i][j] << "    ";
		}
		cout << "\n";
	}

	double res2 = 0.0;

	for(size_t i = 0; i < 3; ++i){
		res2 += A[i][0] * pow(0.333333, 1) + A[i][1] * pow(0.333333, 2) + A[i][2] * pow(0.333333, 3) + A[i][3] * pow(0.333333, 4);
	}

	cout << "\n I = " << res2 << "\n";

	// double res1 = 0.0;
	// for (size_t i = 0; i < n - 1; ++i)
	// {
	// 	if (x[i + 1] - xx > 0 && xx - x[i] > 0)
	// 	{
	// 		res1 = A[i][0] + A[i][1] * (xx - x[i]) + A[i][2] * pow(xx - x[i], 2) + A[i][3] * pow(xx - x[i], 3);
	// 		cout << "Spline (uniform grid): \n";
	// 		cout << res1 << "\n";
	// 		cout << res1 << "\n";
	// 		cout << fabs(res - xx) << "\n";
	// 		// cout << fabs(res - xx*xx) << "\n";
	// 		break;
	// 	}
	// }

	return 0;
}
