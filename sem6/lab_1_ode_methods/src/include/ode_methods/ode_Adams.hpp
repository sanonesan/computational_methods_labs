#pragma once

#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>

#include"vector_norm.hpp"
#include"ode_RK4.hpp"

#define eps 1e-16


/*
	Метод Адамса-Башфорта (первые 3 шага --- Рунге-Кутта 4 порядка)
*/
template<typename T, typename F>
void ode_AB4(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const std::string &out_path){

    std::ofstream fout(out_path);
	if(!fout){
		std::cout << "\n error \n";
		return;
	}

	fout << std::scientific;
	fout << "time"; 
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << ",u" << i;
	}
	fout << "\n";
	fout << start_time;	
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << "," << x[i];
	}
	fout << "\n";


	
	std::vector<T> tmp;
    std::vector<std::vector<T>> func_value;

    for(std::size_t i = 0; i < func.size(); ++i){
        tmp.push_back(func[i](x, start_time));
        func_value.push_back(tmp);
        tmp.clear();
    }
    tmp.assign(x.begin(), x.end());

    //first 3 steps with ode_RK4
    for(std::size_t j = 0; j < 3; ++j){

        fout << start_time + tau;
		for(std::size_t i = 0; i < func.size(); ++i){
			x[i] += runge_coef(start_time, tau, tmp, i, func[i]);
			fout << "," << x[i];
		}
		fout << "\n";

		tmp.assign(x.begin(), x.end());
		start_time += tau;

        for(std::size_t i = 0; i < func.size(); ++i){            
            func_value[i].push_back(func[i](x, start_time));
        }

    }

	while (start_time <= end_time){
		
		start_time += tau;
        
        fout << start_time;
        for(std::size_t i = 0; i < func.size(); ++i){
			x[i] += tau * (2.29167 * func_value[i][3] - 2.45833 * func_value[i][2] + 1.54167 * func_value[i][1] - 0.375 * func_value[i][0]);
            fout << "," << x[i];
		}
        fout << "\n";

        for(std::size_t i = 0; i < func.size(); ++i){
            func_value[i].erase(func_value[i].begin());           
            func_value[i].push_back(func[i](x, start_time));
        }

	}

	fout.close();

    return;
}


/*
	Метод прогноз-коррекция (первые 3 шага --- Рунге-Кутта 4 порядка)
*/
template<typename T, typename F>
void ode_Predictor_Corrector(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const std::string &out_path){

    std::ofstream fout(out_path);
	if(!fout){
		std::cout << "\n error \n";
		return;
	}

	fout << std::scientific;
	fout << "time"; 
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << ",u" << i;
	}
	fout << "\n";
	fout << start_time;
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << "," << x[i];
	}
	fout << "\n";
	
	std::vector<T> tmp;
    std::vector<std::vector<T>> func_value;

    for(std::size_t i = 0; i < func.size(); ++i){
        tmp.push_back(func[i](x, start_time));
        func_value.push_back(tmp);
        tmp.clear();
    }
    tmp.assign(x.begin(), x.end());

    //first 3 steps with ode_RK4
    for(std::size_t j = 0; j < 3; ++j){

        fout << start_time + tau;
		for(std::size_t i = 0; i < func.size(); ++i){
			x[i] += runge_coef(start_time, tau, tmp, i, func[i]);
			fout << "," << x[i];
		}
		fout << "\n";

		tmp.assign(x.begin(), x.end());
		start_time += tau;

        for(std::size_t i = 0; i < func.size(); ++i){            
            func_value[i].push_back(func[i](x, start_time));
        }

    }

	while (start_time <= end_time){
		
        //prediction (AB4 step)
		start_time += tau;
        for(std::size_t i = 0; i < func.size(); ++i){
			x[i] += tau * (2.29167 * func_value[i][3] - 2.45833 * func_value[i][2] + 1.54167 * func_value[i][1] - 0.375 * func_value[i][0]);
		}

        for(std::size_t i = 0; i < func.size(); ++i){
            func_value[i].erase(func_value[i].begin());           
            func_value[i].push_back(func[i](x, start_time));
        }

        //correction
        fout << start_time;
        for(std::size_t i = 0; i < func.size(); ++i){
			x[i] += tau * (0.375 * func_value[i][3] + 0.791667 * func_value[i][2] - 0.208333 * func_value[i][1] + 0.0416667 * func_value[i][0]);
            fout << "," << x[i];
        }
        fout << "\n";

	}

	fout.close();

    return;
}
