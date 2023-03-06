#pragma once

#include<vector>
#include<cmath>
#include <functional> 


template<typename T>
class Var_4{

	private:

    public:

	std::vector<T> _x0 = {1.,0.1};
	std::vector<std::function<T (const std::vector<T>& x, const T t)>> _ode_system;


    Var_4(){
	

		const T a = 0.1, mu = 0.1, w = 0.25, nu = 0.15;

		//H
		auto func_H = [a](const std::vector<T> &x) -> T{
			return pow(x[0] * x[0] + x[1] * x[1], 2) - 2 * a * a * (x[0] * x[0] - x[1] * x[1]);
		};


		//num_deriv_dH_dx_i
		auto func_num_deriv_dH_dx_i = [func_H](std::vector<T> x, const size_t i) -> T{
			
			if(i >= x.size())
				throw std::invalid_argument("vector out of range");
			T res = - func_H(x);
			T _eps = 1e-9;
			x[i] += _eps;
			res += func_H(x);
			return res / _eps;
		};


		//exact_deriv_dH_dx_i
		auto func_exact_deriv_dH_dx_i = [a, func_H](std::vector<T> x, const std::size_t i) -> T{
			
			if(i == 0)
				return -4*a*a*x[0] + 4 * x[0] * (x[0]*x[0] + x[1]*x[1]);
			if(i == 1)
					return 4*a*a*x[1] + 4 * x[1] * (x[0]*x[0] + x[1]*x[1]);
			else{
				throw std::invalid_argument("vector out of range");
			}
		};


		// // dx/dt
		// auto func_1 = [mu, nu, w, func_H, func_num_deriv_dH_dx_i](const std::vector<T>& x, const T t) -> T{
		// 	return func_num_deriv_dH_dx_i(x, 1) - mu * func_H(x) * func_num_deriv_dH_dx_i(x, 0) + nu * x[1] * sin(w * t);
		// };

		// //dy/dt
		// auto func_2 = [mu, nu, w, func_H, func_num_deriv_dH_dx_i](const std::vector<T>& x, const T t) -> T{
		// 	return -func_num_deriv_dH_dx_i(x, 0) - mu * func_H(x) * func_num_deriv_dH_dx_i(x, 1) + nu * x[1] * sin(w * t);
		// };


		// dx/dt
		auto func_1 = [mu, nu, w, func_H, func_exact_deriv_dH_dx_i](const std::vector<T>& x, const T t) -> T{
			return func_exact_deriv_dH_dx_i(x, 1) - mu * func_H(x) * func_exact_deriv_dH_dx_i(x, 0) + nu * x[1] * sin(w * t);
		};


		// dy/dt
		auto func_2 = [mu, nu, w, func_H, func_exact_deriv_dH_dx_i](const std::vector<T>& x, const T t) -> T{
			return -func_exact_deriv_dH_dx_i(x, 0) - mu * func_H(x) * func_exact_deriv_dH_dx_i(x, 1) + nu * x[1] * sin(w * t);
		};


		_ode_system.push_back(func_1);
		_ode_system.push_back(func_2);

		}

};
