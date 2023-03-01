#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include <functional>
#define eps 1e-9

using namespace std;

template<typename T>
void print_vec(const vector<T>& vec){
	cout << "\n( ";
	for(size_t i = 0; i < vec.size() - 1; ++i)
		cout << vec[i] << "\t";
	cout << vec[vec.size()-1] << " )^T \n";
}

template<typename T>
T norm(vector<T> vec){
	return sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
}

double a = 0.1, mu = 0.1, w = 0.25, nu = 0.15;

template<typename T, typename F>
T runge_coef(const T t, const T tau, const vector<T>& x, const size_t k, const F &func){

	vector<T> tmp(x);
	T k1, k2, k3, k4;

	k1 = tau * func(tmp, t);
	tmp[k] = x[k] + 0.5 * k1;

	k2 = tau * func(tmp, t + 0.5 * tau);
	tmp[k] = x[k] + 0.5 * k2;

	k3 = tau * func(tmp, t + 0.5 * tau);
	tmp[k] = x[k] + k3;

	k4 = tau * func(tmp, t + tau);

	return 0.166666 * (k1 + k4) + 0.333333 * (k2 + k3);
}

template<typename T, typename F>
void runge_cutta(T start_time, T end_time, T tau, vector<T> x, const vector<F> &func, const string &out_path){

	ofstream fout(out_path);
	if(!fout){
		cout << "\n error \n";
	}

	vector<T> tmp(x);
	T coef = 0.;

	while (start_time <= end_time){

		for(size_t i = 0; i < func.size(); ++i){
			x[i] += runge_coef(start_time, tau, tmp, i, func[i]);
			fout << "\t" << x[i];
		}
		for(size_t i = 0; i < func.size(); ++i)
			tmp[i] = x[i];

		fout << "\n";

		start_time += tau;
	}

	fout.close();
}

template<typename T, typename F>
void runge_cutta_vary_step(T start_time, T end_time, T tau, vector<T> x, const vector<F> &func, const string &out_path){

	ofstream fout(out_path);
	if(!fout){
		cout << "\n error \n";
	}

	vector<T> tmp(x);
	vector<T> x_1(x);
	vector<T> x_2(x);
	bool flag = true;
	while (start_time <= end_time){

		for(size_t i = 0; i < func.size(); ++i){
			x_1[i] += runge_coef(start_time, tau, x, i, func[i]);
			x_2[i] += runge_coef(start_time, tau / 2, x, i, func[i]);
			x_2[i] += runge_coef(start_time, tau / 2, x_2, i, func[i]);
			
		}
		for(size_t i = 0; i < func.size(); ++i){
			tmp[i] = x_1[i] - x_2[i];
		}
		while (norm(tmp) / 15 >= 0.001 && tau > 0.002)
		{
			tau /= 2;
			for(size_t i = 0; i < func.size(); ++i){
				x_1[i] += runge_coef(start_time, tau, x, i, func[i]);
				x_2[i] += runge_coef(start_time, tau / 2, x, i, func[i]);
				x_2[i] += runge_coef(start_time, tau / 2, x_2, i, func[i]);
			}
			for(size_t i = 0; i < func.size(); ++i){
				tmp[i] = x_1[i] - x_2[i];
		}
		}	

		if (norm(tmp) / 15 < 0.001 || tau < 0.002)
			tau *= 2;

		for(size_t i = 0; i < func.size(); ++i){
			x[i] = x_1[i];
			fout << "\t" << x_1[i];
		}
		


		fout << "\n";
		start_time += tau;
	}

	fout.close();
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

	//exact_deriv_dH_dx_i
	auto func_dH_dx_i_exact = [func_H](vector<T> x, const size_t i) -> T{
		if(i == 0)
			return -4*a*a*x[0] + 4 * x[0] * (x[0]*x[0] + x[1]*x[1]);
		if(i == 1)
				return 4*a*a*x[1] + 4 * x[1] * (x[0]*x[0] + x[1]*x[1]);
		else{
			throw __throw_logic_error;
		}
	};

	// // dx/dt
	// auto func_1 = [func_H, func_dH_dx_i](const vector<T>& x, const T t) -> T{
	// 	return func_dH_dx_i(x, 1) - mu * func_H(x) * func_dH_dx_i(x, 0) + nu * x[1] * sin(w * t);
	// };

	// //dy/dt
	// auto func_2 = [func_H, func_dH_dx_i](const vector<T>& x, const T t) -> T{
	// 	return -func_dH_dx_i(x, 0) - mu * func_H(x) * func_dH_dx_i(x, 1) + nu * x[1] * sin(w * t);
	// };


	// dx/dt
	auto func_1 = [func_H, func_dH_dx_i_exact](const vector<T>& x, const T t) -> T{
		return func_dH_dx_i_exact(x, 1) - mu * func_H(x) * func_dH_dx_i_exact(x, 0) + nu * x[1] * sin(w * t);
	};

	// dy/dt
	auto func_2 = [func_H, func_dH_dx_i_exact](const vector<T>& x, const T t) -> T{
		return -func_dH_dx_i_exact(x, 0) - mu * func_H(x) * func_dH_dx_i_exact(x, 1) + nu * x[1] * sin(w * t);
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

	t = 0;
	t_final = 50;
	tau = 0.01;
	x[0] = 1.0;
	x[1] = 0.1;

	out_path = "RK4_vary_output.txt";
	runge_cutta_vary_step(t, t_final, tau, x, _functions, out_path);


    return 0;
}

