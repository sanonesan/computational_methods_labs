#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include <functional>
#define eps 1e-9

using namespace std;

double a = 0.1, mu = 0.1, w = 0.25, nu = 0.15;

template<typename T, typename F>
T runge_coef(const T t, const T tau, const vector<T>& x, F func){

	vector<T> tmp(x);

	T k1, k2, k3, k4;

	k1 = tau * func(tmp, t);
	for(size_t k = 0; k < x.size(); ++k)
		tmp[k] = x[k] + 0.5 * k1;	
	k2 = tau * func(tmp, t + 0.5 * tau);
	for(size_t k = 0; k < x.size(); ++k)
		tmp[k] = x[k] + 0.5 * k2;
	k3 = tau * func(tmp, t + 0.5 * tau);
	for(size_t k = 0; k < x.size(); ++k)
		tmp[k] = x[k] + k3;
	k4 = tau * func(tmp, t + tau);
	
	return 0.166666 * (k1 + k4) + 0.333333 * (k2 + k3);
}

template<typename T, typename F>
T runge_cutta(T start_time, T end_time, T tau, vector<T> x, vector<F> func, string out_path){

	ofstream fout(out_path);
	if(!fout){
		cout << "\n error \n";
	}

	vector<T> tmp(x);

	//fout << scientific;
	//fout.precision(8);

	//fout << "%time(s) \t x-pos \t y-pos \n" ;
	//fout << start_time << "\t" << x[0] << "\t"<<  x[1] << "\n";
	
	T coef = 0.;

	while (start_time <= end_time){

		//fout << start_time; 
		for(size_t i = 0; i < func.size(); ++i){
			x[i] += runge_coef(start_time, tau, tmp, func[i]);
			fout << "\t" << x[i];
		}
		for(size_t i = 0; i < func.size(); ++i)
			tmp[i] = x[i];
		
		fout << "\n";	

		start_time += tau;
	}

	fout.close();
		
}

template<typename T>
void print_vec(const vector<T>& vec){
	cout << "\n( ";
	for(size_t i = 0; i < vec.size() - 1; ++i)
		cout << vec[i] << "\t";
	cout << vec[vec.size()-1] << " )^T \n";
}



int main(int args, char **argv){

	typedef double T;

    T t, t_final, tau;

	
	// ---------------------------- //
	// ---------funcs_var_4-------- //
	// ---------------------------- //

	//H
	auto func_H = [](const vector<T> &x) -> T{
		return pow(x[0] * x[0] + x[1] * x[1], 2) - 2 * a * a * (x[0] * x[0] - x[1] * x[1]);
	};

	//num_deriv_dH_dx_i
	auto func_dH_dx_i = [func_H](vector<T> x, const size_t i) -> T{
		T res = - func_H(x);
		x[i] += eps;
		res += func_H(x);
		return res / eps;
	};

	// dx/dt
	auto func_1 = [func_H, func_dH_dx_i](const vector<T>& x, const T t) -> T{
		return func_dH_dx_i(x, 1) - mu * func_H(x) * func_dH_dx_i(x, 0) + nu * x[1] * sin(w * t);
	};

	// dy/dt
	auto func_2 = [func_H, func_dH_dx_i](const vector<T>& x, const T t) -> T{
		return -func_dH_dx_i(x, 0) - mu * func_H(x) * func_dH_dx_i(x, 1) + nu * x[1] * sin(w * t);
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

	string out_path = "RK4_output.txt";

	runge_cutta(t, t_final, tau, x, _functions, out_path);


    return 0;
}

