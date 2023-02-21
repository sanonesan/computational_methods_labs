#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#define eps 1e-9

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;

double a = 0.1, mu = 0.1, w = 0.25, nu = 0.15;

//H
template<typename T>
T func_H(const T x, const T y){
	return pow(x * x + y * y, 2) - 2 * a * a * (x * x - y * y);
}

//num_deriv_dH_dx
template<typename T>
T func_dH_dx(const T x, const T y){

	return (func_H(x + eps, y) - func_H(x, y)) / eps;
}

//num_deriv_dH_dy
template<typename T>
T func_dH_dy(const T x, const T y){
	return (func_H(x , y + eps) - func_H(x, y)) / eps;
}

// dx/dt
template<typename T>
T func_1(const T x, const T y, const T t){
	return func_dH_dy(x, y) - mu * func_H(x, y) * func_dH_dx(x, y) + nu * y * sin(w * t);
}

// dy/dt
template<typename T>
T func_2(const T x, const T y, const T t){
	return -func_dH_dx(x, y) - mu * func_H(x, y) * func_dH_dy(x, y) + nu * y * sin(w * t);
}




int main(int args, char **argv){

    typedef double Type;

    Type t, t_final, tau;

	vector<Type> x,y;

	t = 0;
	t_final = 50;
	tau = 0.01;

	x.push_back(1.0);
	y.push_back(0.1);

	ofstream fout("RK4_output.csv");
	if(!fout){
		cout << "\n error \n";
	}

	fout << scientific;
	fout.precision(8);

	fout << "%time(s) \t x-pos \t y-pos \n" ;

	Type k1, k2, k3, k4;
	Type l1, l2, l3, l4;
	Type t1, t2, t3, t4;
	int counter = 0;
	while (t <= t_final){
		
		k1 = tau * func_1(x[counter], y[counter], t);
		l1 = tau * func_2(x[counter], y[counter], t);

		k2 = tau * func_1(x[counter] + 0.5 * k1, y[counter] + 0.5 * l1, t);
		l2 = tau * func_2(x[counter] + 0.5 * k1, y[counter] + 0.5 * l1, t);

		k3 = tau * func_1(x[counter] + 0.5 * k2, y[counter] + 0.5 * l2, t);
		l3 = tau * func_2(x[counter] + 0.5 * k2, y[counter] + 0.5 * l2, t);

		k4 = tau * func_1(x[counter] + 0.5 * k3, y[counter] + 0.5 * l3, t);
		l4 = tau * func_2(x[counter] + 0.5 * k3, y[counter] + 0.5 * l3, t);

		x.push_back(x[counter] + 0.166 * k1 + 0.333 * k2 + 0.333 * k3 + 0.166 * k4);
		y.push_back(y[counter] + 0.166 * l1 + 0.333 * l2 + 0.333 * l3 + 0.166 * l4);


		++counter;
		t+=tau;

		fout << t << "\t" << x[counter] << "\t"<<  y[counter] << "\n";

	}

    return 0;
}

