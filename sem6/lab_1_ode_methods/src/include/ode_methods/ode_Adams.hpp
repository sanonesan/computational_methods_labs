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
	fout << std::setprecision(48);
	fout << "time"; 
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << ",u" << i;
	}
	fout << ",tau\n";
	fout << start_time;	
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << "," << x[i];
	}
	fout << "," << tau << "\n";



	
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
			x[i] += RK4_coef(start_time, tau, tmp, i, func[i]);
			fout << "," << x[i];
		}
		fout << "," << tau << "\n";


		tmp.assign(x.begin(), x.end());
		start_time += tau;

        for(std::size_t i = 0; i < func.size(); ++i){            
            func_value[i].push_back(func[i](x, start_time));
        }

    }

	T c1 = 55;
	T c2 = 59;
	T c3 = 37;
	T c4 = 0.375;

	c1 /= 24;
	c2 /= 24;
	c3 /= 24;
	
	while (start_time <= end_time){
		
		start_time += tau;
        
        fout << start_time;
        for(std::size_t i = 0; i < func.size(); ++i){
			x[i] += tau * (c1 * func_value[i][3] - c2 * func_value[i][2] + c3 * func_value[i][1] - c4 * func_value[i][0]);
            fout << "," << x[i];
		}
		fout << "," << tau << "\n";


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
	fout << std::setprecision(48);
	fout << "time"; 
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << ",u" << i;
	}
	fout << ",tau\n";
	fout << start_time;
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << "," << x[i];
	}
	fout << "," << tau << "\n";
	
	
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
			x[i] += RK4_coef(start_time, tau, tmp, i, func[i]);
			fout << "," << x[i];
		}
		fout << "," << tau << "\n";

		tmp.assign(x.begin(), x.end());
		start_time += tau;

        for(std::size_t i = 0; i < func.size(); ++i){            
            func_value[i].push_back(func[i](x, start_time));
        }

    }

	T c1 = 55;
	T c2 = 59;
	T c3 = 37;
	T c4 = 0.375;
	T c5 = 19;
	T c6 = 5;

	c1 /= 24;
	c2 /= 24;
	c3 /= 24;
	c5 /= 24;
	c6 /= 24;



	while (start_time <= end_time){

		tmp.assign(x.begin(), x.end());
		
        //prediction (AB4 step)
		start_time += tau;
        for(std::size_t i = 0; i < func.size(); ++i){
			tmp[i]  += tau * (c1 * func_value[i][3] - c2 * func_value[i][2] + c3 * func_value[i][1] - c4 * func_value[i][0]);
		}

        for(std::size_t i = 0; i < func.size(); ++i){
            func_value[i].erase(func_value[i].begin());           
            func_value[i].push_back(func[i](tmp, start_time));
        }

        //correction
        fout << start_time;
        for(std::size_t i = 0; i < func.size(); ++i){
			x[i] += tau * (c4 * func_value[i][3] + c5 * func_value[i][2] - c6 * func_value[i][1] + func_value[i][0] / 24);
            fout << "," << x[i];
        }
		fout << "," << tau << "\n";

		for(std::size_t i = 0; i < func.size(); ++i){
            func_value[i][3] = func[i](x, start_time);
        }
	}

	fout.close();

    return;
}
