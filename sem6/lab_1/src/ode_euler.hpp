#pragma once

#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>

#include"vector_norm.hpp"

#define eps 1e-16

using namespace std;


/*
    Явный метод Эйлера
*/
template<typename T, typename F>
void ode_euler(T start_time, T end_time, T tau, vector<T> x, const vector<F> &func, const string &out_path){

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

    while(start_time <= end_time){  
        fout << start_time + tau;
        for(size_t i = 0; i < x.size(); ++i){
            x[i] = tmp[i] + tau * func[i](tmp, start_time);
            fout << "," << x[i];
        }
        fout << "\n";
        tmp.assign(x.begin(), x.end());
        start_time += tau;
    }    

    return;
}