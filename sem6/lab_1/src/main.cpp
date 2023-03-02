#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include <functional>

using namespace std;

#include"RK4.hpp"


int main(int args, char **argv){

	typedef double T;

	T tol = 0.001; //tolerance
    T t, t_final, tau;

	// ---------------------------- //
	// ---------funcs_var_4-------- //
	// ---------------------------- //

	//const
	const T a = 0.1, mu = 0.1, w = 0.25, nu = 0.15;

	//H
	auto func_H = [a](const vector<T> &x) -> T{
		return pow(x[0] * x[0] + x[1] * x[1], 2) - 2 * a * a * (x[0] * x[0] - x[1] * x[1]);
	};

	//num_deriv_dH_dx_i
	auto func_num_deriv_dH_dx_i = [func_H](vector<T> x, const size_t i) -> T{
		
		if(i >= x.size())
			throw invalid_argument("vector out of range");
		T res = - func_H(x);
		T _eps = 1e-9;
		x[i] += _eps;
		res += func_H(x);
		return res / _eps;
	};

	//exact_deriv_dH_dx_i
	auto func_exact_deriv_dH_dx_i = [a, func_H](vector<T> x, const size_t i) -> T{
		
		if(i == 0)
			return -4*a*a*x[0] + 4 * x[0] * (x[0]*x[0] + x[1]*x[1]);
		if(i == 1)
				return 4*a*a*x[1] + 4 * x[1] * (x[0]*x[0] + x[1]*x[1]);
		else{
			throw invalid_argument("vector out of range");
		}
	};

	// // dx/dt
	// auto func_1 = [mu, nu, w, func_H, func_num_deriv_dH_dx_i](const vector<T>& x, const T t) -> T{
	// 	return func_num_deriv_dH_dx_i(x, 1) - mu * func_H(x) * func_num_deriv_dH_dx_i(x, 0) + nu * x[1] * sin(w * t);
	// };

	// //dy/dt
	// auto func_2 = [mu, nu, w, func_H, func_num_deriv_dH_dx_i](const vector<T>& x, const T t) -> T{
	// 	return -func_num_deriv_dH_dx_i(x, 0) - mu * func_H(x) * func_num_deriv_dH_dx_i(x, 1) + nu * x[1] * sin(w * t);
	// };


	// dx/dt
	auto func_1 = [mu, nu, w, func_H, func_exact_deriv_dH_dx_i](const vector<T>& x, const T t) -> T{
		return func_exact_deriv_dH_dx_i(x, 1) - mu * func_H(x) * func_exact_deriv_dH_dx_i(x, 0) + nu * x[1] * sin(w * t);
	};

	// dy/dt
	auto func_2 = [mu, nu, w, func_H, func_exact_deriv_dH_dx_i](const vector<T>& x, const T t) -> T{
		return -func_exact_deriv_dH_dx_i(x, 0) - mu * func_H(x) * func_exact_deriv_dH_dx_i(x, 1) + nu * x[1] * sin(w * t);
	};

	vector<function<T (const vector<T>& x, const T t)>> _functions;
	vector<T> x = {1.,0.1};
	_functions.push_back(func_1);
	_functions.push_back(func_2);


	// ---------------------------- //
	// ---------funcs_var_4-------- //
	// ---------------------------- //


	t = 0;
	t_final = 50;
	tau = 0.01;

	// string out_path = "../output/RK4_fix_step_output.txt";
	string out_path = "../output/RK4_fix_step_output.csv";
	runge_cutta_fix_step(t, t_final, tau, x, _functions, out_path);

	tol = 0.001;
	t = 0;
	t_final = 50;
	tau = 0.001;
	x[0] = 1.0;
	x[1] = 0.01;

	out_path = "../output/RK4_vary_step_output.csv";
	runge_cutta_vary_step(t, t_final, tau, x, _functions, tol, out_path);

    return 0;
}

