#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>

#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "../Class_1d_wave_equation.hpp"

template <typename T>
void explicit_cross_difference_scheme(

    const Class_1d_wave_equation<T>& wave_equation,
    const std::vector<T>& time,
    const std::vector<T>& x,
    std::vector<T>& y_j_minus_1,
    const T& tolerance,
    std::ofstream &fout){
    

    std::vector<T> y_j(y_j_minus_1);
    std::vector<T> y_j_plus_1(y_j_minus_1);

    const T a_tau_h_square = pow(wave_equation._a, 2) * pow(wave_equation._tau, 2) / pow(wave_equation._h, 2);
    const T a_tau_h_square__2 = a_tau_h_square / 2;

    y_j[0] = y_j_minus_1[0];
    for(std::size_t i = 1; i < x.size() - 1; ++i){
        y_j[i] = y_j_minus_1[i];
        y_j[i] += wave_equation._tau * wave_equation._initial_conditions[1](x[i], time[0]);
        y_j[i] += a_tau_h_square__2 * (y_j_minus_1[i+1] - 2 * y_j_minus_1[i] + y_j_minus_1[i-1]);
    }
    y_j[y_j.size() - 1] = y_j_minus_1[y_j.size() - 1];

    // y_j[0] = y_j_minus_1[0];
    // for(std::size_t i = 1; i < x.size() - 1; ++i){
    //     y_j[i] = y_j_minus_1[i];
    //     y_j[i] += wave_equation._tau * wave_equation._initial_conditions[1](x[i], time[0]);
    //     y_j[i] += a_tau_h_square__2 * wave_equation._initial_conditions[2](x[i], time[0]);//(y_j_minus_1[i+1] - 2 * y_j_minus_1[i] + y_j_minus_1[i-1]);
    // }
    // y_j[y_j.size() - 1] = y_j_minus_1[y_j.size() - 1];

    fout << y_j[0];
    for (std::size_t i = 1; i < x.size() - 1; ++i) {
        fout << "," << y_j[i];
    }
    fout << "," << y_j[y_j.size() - 1] << "\n";
    
    for (std::size_t j = 2; j < time.size(); ++j) {

        y_j_plus_1[0] = wave_equation._boundary_conditions[0](x[0], time[j]);
        for (std::size_t i = 1; i < y_j.size() - 1; ++i) {
            y_j_plus_1[i] = a_tau_h_square * (y_j[i + 1] - 2 * y_j[i] + y_j[i - 1]) + 2 * y_j[i] - y_j_minus_1[i];
        }
        y_j_plus_1[y_j.size() - 1] = wave_equation._boundary_conditions[1](x[x.size() - 1], time[j]);

        y_j_minus_1.assign(y_j.begin(), y_j.end());
        y_j.assign(y_j_plus_1.begin(), y_j_plus_1.end());

        fout << y_j[0];
        for (std::size_t i = 1; i < x.size() - 1; ++i) {
            fout << "," << y_j[i];
        }
        fout << "," << y_j[y_j.size() - 1] << "\n";
    }
};
