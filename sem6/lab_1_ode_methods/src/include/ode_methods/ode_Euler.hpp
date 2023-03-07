#pragma once

#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>

#define eps 1e-16


/*
    Явный метод Эйлера
*/
template<typename T, typename F>
void ode_Euler(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const std::string &out_path){

    std::ofstream fout(out_path);
	if(!fout){
		std::cout << "\n error \n";
		return;
	}

	fout << std::scientific;
	fout << "time"; 
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << "," << "u" << i;
	}
	fout << "\n";
	fout << start_time;
	for(std::size_t i = 0; i < func.size(); ++i){
		fout << "," << x[i];
	}
	fout << "\n";
	
	std::vector<T> tmp(x);

    while(start_time <= end_time){  
        fout << start_time + tau;
        for(std::size_t i = 0; i < x.size(); ++i){
            x[i] = tmp[i] + tau * func[i](tmp, start_time);
            fout << "," << x[i];
        }
        fout << "\n";
        tmp.assign(x.begin(), x.end());
        start_time += tau;
    }    

    return;
}