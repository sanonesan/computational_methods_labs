#pragma once

#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <vector>
#include "../Class_1d_heat_equation.hpp"


template <typename T>
void scheme_1(Class_1d_heat_equation<T>& heat_equation, const std::string &out_path) {
    
    std::ofstream fout(out_path + "_t.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
	fout << std::setprecision(8);
    
    std::vector<T> time;
    time.push_back(heat_equation._start_time);
    fout << time[time.size() - 1] << "\n";
    while(time[time.size() - 1] <= heat_equation._end_time){
        time.push_back(time[time.size() - 1] + heat_equation._tau);
        fout << time[time.size() - 1] << "\n";
    }

    fout.close();

    fout.open(out_path + "_x.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
	fout << std::setprecision(8);

    std::vector<T> x;
    x.push_back(heat_equation._x0);
    fout << x[x.size() - 1] << "\n";
    while(x[x.size() - 1] + heat_equation._h <= heat_equation._xL){
        x.push_back(x[x.size() - 1] + heat_equation._h);
        fout << x[x.size() - 1] << "\n";
    }

    fout.close();

    fout.open(out_path + "_y.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
	fout << std::setprecision(8);

    std::vector<T> y;
    y.assign(x.begin(), x.end());

    y[0] = heat_equation._boundary_conditions[0](x[0], time[0]);
    fout << y[0];
    for(std::size_t i = 1; i < x.size()-1; ++i){
        y[i] = heat_equation._initial_conditions(x[i], time[0]);
        fout << "," << y[i];

    }
    y[y.size() - 1] = heat_equation._boundary_conditions[1](x[x.size() - 1], time[0]);
    fout << "," << y[y.size() - 1] << "\n";
    
    std::vector<T> y_prev;
    y_prev.assign(y.begin(), y.end());

    std::vector<T> a_i;
    a_i.assign(x.begin(), x.end());
    for(std::size_t i = 1; i < x.size(); ++i){
        a_i[i] = heat_equation._K(x[i] - heat_equation._h / 2);
    }

    T C = 0.;
    C = heat_equation._tau / ( pow(heat_equation._h, 2) * heat_equation._c * heat_equation._rho );

    for(std::size_t j = 1; j < time.size(); ++j){

        y[0] = heat_equation._boundary_conditions[0](x[0], time[j]);
        fout << y[0];
        for(std::size_t i = 1; i < y.size() - 1; ++i){
            y[i] = C * ( a_i[i + 1] * (y_prev[i+1] - y_prev[i]) - a_i[i] * (y_prev[i] - y_prev[i-1])) + y_prev[i];
            fout << "," << y[i];
        }
        y[y.size() - 1] = heat_equation._boundary_conditions[1](x[x.size() - 1], time[j]);
        fout << "," << y[y.size() - 1] << "\n";
        y_prev.assign(y.begin(), y.end());

    }


    fout.close();

};
