#pragma once

#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>

#include"vector_norm.hpp"

#define eps 1e-16

using namespace std;


/* 
    template<typename T>
    T runge_coef(const T t, const T tau, const vector<T>& x, const size_t k, const F &func);

    template<typename T, typename F>
    void RK4_fix_step(T start_time, T end_time, T tau, vector<T> x, const vector<F> &func, const string &out_path);

    template<typename T, typename F>
    void RK4_vary_step(T start_time, T end_time, T tau, vector<T> x, const vector<F> &func, const T tol, const string &out_path);
*/


/*
	Вычисление коэффициентов метода Рунге-Кутты
*/
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

/*
	Метод Рунге-Кутты 4 порядка с фиксированным шагом
*/
template<typename T, typename F>
void RK4_fix_step(T start_time, T end_time, T tau, vector<T> x, const vector<F> &func, const string &out_path){

	ofstream fout(out_path);
	if(!fout){
		cout << "\n error \n";
		return;
	}

	fout << scientific;
	fout << "time," << "x," << "y\n";
	fout << start_time;
	for(size_t i = 0; i < func.size(); ++i){
		fout << "," << x[i];
	}
	fout << "\n";
	
	vector<T> tmp(x);

	while (start_time <= end_time){
		
		fout << start_time + tau;

		for(size_t i = 0; i < func.size(); ++i){
			x[i] += runge_coef(start_time, tau, tmp, i, func[i]);
			fout << "," << x[i];
		}
		fout << "\n";
		
		tmp.assign(x.begin(), x.end());

		start_time += tau;
	}

	fout.close();

    return;
}

/*
	Метод Рунге-Кутты 4 порядка с изменяющимя шагом
*/
template<typename T, typename F>
void RK4_vary_step(T start_time, T end_time, T tau, vector<T> x, const vector<F> &func, const T tol, const string &out_path){

	ofstream fout(out_path);
	if(!fout){
		cout << "\n error \n";
		return;
	}

	fout << scientific;
	fout << "time," << "x," << "y\n";
	fout << start_time;
	for(size_t i = 0; i < func.size(); ++i){
		fout << "," << x[i];
	}
	fout << "\n";


	vector<T> tmp(x);
	vector<T> x_1(x);
	vector<T> x_2(x);
	T breaker;
	T coef_tol = tol / 5 / 10000;

	while (true){

		if (start_time + tau > end_time){
			
			if((start_time < end_time)){

				tau = end_time - start_time;
				fout << start_time + tau;
				for(size_t i = 0; i < func.size(); ++i){
					x[i] += runge_coef(start_time, tau, x_1, i, func[i]);
					fout << "," << x[i];
				}
				fout << "\n";
			}
			break;
		}

		while (true) {
			
			for(size_t i = 0; i < func.size(); ++i){
				x_1[i] = x[i] + runge_coef(start_time, tau, x, i, func[i]);
				tmp[i] = x[i] + runge_coef(start_time, tau / 2, x, i, func[i]);
			}
			for(size_t i = 0; i < func.size(); ++i)
				x_2[i] = tmp[i] + runge_coef(start_time, tau / 2, tmp, i, func[i]);

			for(size_t i = 0; i < func.size(); ++i){
				tmp[i] = x_1[i] - x_2[i];
			}		

			breaker = norm_inf(tmp) / 15;
			if(breaker >= tol){
				tau /= 2;
			}
			else{
				x.assign(x_1.begin(), x_1.end());
				break;
			}

		}		

		start_time += tau;	

		if (breaker < coef_tol)
			tau *= 2;

		fout << start_time;
		for(size_t i = 0; i < func.size(); ++i){
			fout << "," << x[i];
		}
		fout << "\n";

	}

	fout.close();

    return;
}
