#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "../Class_1d_heat_equation.hpp"
#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"


template <typename T>
void implicit_2_layer_difference_scheme(const Class_1d_heat_equation<T>& heat_equation, const std::string &out_path) {
    
    //time steps + output

    std::ofstream fout(out_path + "_t.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
	fout << std::setprecision(8);

    //init time vector
    std::vector<T> time;

    fout << "time\n";
    time.push_back(heat_equation._start_time);
    fout << time[time.size() - 1] << "\n";
    while(time[time.size() - 1] <= heat_equation._end_time){
        time.push_back(time[time.size() - 1] + heat_equation._tau);
        fout << time[time.size() - 1] << "\n";
    }

    fout.close();

    //space steps + output

    fout.open(out_path + "_x.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
	fout << std::setprecision(8);

    //init space vector
    std::vector<T> x;

    fout << "x\n";
    x.push_back(heat_equation._x0);
    fout << x[x.size() - 1] << "\n";
    while(x[x.size() - 1] + heat_equation._h <= heat_equation._xL){
        x.push_back(x[x.size() - 1] + heat_equation._h);
        fout << x[x.size() - 1] << "\n";
    }

    fout.close();


    //Solution

    fout.open(out_path + "_y.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
	fout << std::setprecision(8);

    //init solution vector
    std::vector<T> y;
    y.assign(x.begin(), x.end());
    
    //output
    fout << "t0";
    for(std::size_t i = 1; i < x.size(); ++i){
        fout << ",t" << i;
    }
    fout << "\n";

    //get initial values for solution vector + output
    y[0] = heat_equation._boundary_conditions[0](x[0], time[0]);
    fout << y[0];
    for(std::size_t i = 1; i < x.size()-1; ++i){
        y[i] = heat_equation._initial_conditions(x[i], time[0]);
        fout << "," << y[i];
    }
    y[y.size() - 1] = heat_equation._boundary_conditions[1](x[x.size() - 1], time[0]);
    fout << "," << y[y.size() - 1] << "\n";
    
    //vector for solutions from previous time step
    std::vector<T> y_prev;
    y_prev.assign(y.begin(), y.end());

    //vector for K(u, x)
    std::vector<T> a_i;
    a_i.assign(x.begin(), x.end());
    for(std::size_t i = 1; i < x.size(); ++i){
        a_i[i] = (heat_equation._K(x[i]) + heat_equation._K(x[i-1])) / 2;
    }

    T C = 0.;
    C = ( pow(heat_equation._h, 2) * heat_equation._c * heat_equation._rho ) / heat_equation._tau ;
    
    Solver_SLE<T> solver_SLE;
    Vector<T> solution(y.size());
    Vector<T> b(y.size());
    Matrix<T> banded_matrix(3, y.size());    
    
    //computating solutions + output
    for(std::size_t j = 1; j < time.size(); ++j){

        b.assign(y.begin(), y.end());
        //std::cout << b << "\n";
        b *= C;

        for(std::size_t i = 0; i < y.size(); ++i){
            // ci
            if (i < y.size() - 1)
                banded_matrix[0][i] = - a_i[i + 1];
            // bi
            banded_matrix[1][i] = a_i[i+1] + a_i[i] + C;
            // ai
            if (i > 0)
                banded_matrix[2][i] = - a_i[i];
            
        }

        solution = solver_SLE.solve_banded(1, 1, banded_matrix, b);        

        y[0] = heat_equation._boundary_conditions[0](x[0], time[j]);
        for(std::size_t i = 1; i < y.size()-1; ++i){
            y[i] = solution[i];
        }
        y[y.size() - 1] = heat_equation._boundary_conditions[0](x[x.size() - 1], time[j]);

        fout << y[0];
        for(std::size_t i = 1; i < x.size()-1; ++i){
            fout << "," << y[i];
        }
        fout << "," << y[y.size() - 1] << "\n";

    }


    fout.close();

};
