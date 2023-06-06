#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "../Class_2d_Poisson_equation_in_rectangle.hpp"
#include "./poisson_FDM_scheme.hpp"

template <typename T>
void Poisson_FDM_in_rectangle(
    const Class_2d_Poisson_equation_in_rectangle<T>& poisson_eq,
    const T& tolerance,
    const std::string& out_path) {
    /**
     *
     *
     * GETTING INITIAL VALUES
     *
     *
     */

    // time steps + output

    std::ofstream fout;
    fout << std::scientific;
    fout << std::setprecision(8);

    fout.open(out_path + "_t.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    // init time vector
    std::vector<T> time;

    fout << "time\n";
    time.resize(std::size_t((poisson_eq._end_time - poisson_eq._start_time) / poisson_eq._tau + 1));
    time[0] = poisson_eq._start_time;
    fout << time[0] << "\n";
    for (std::size_t i = 1; i < time.size() - 1; ++i) {
        time[i] = time[i - 1] + poisson_eq._tau;
        fout << time[i] << "\n";
    }
    time[time.size() - 1] = poisson_eq._end_time;
    fout << time[time.size() - 1] << "\n";

    fout.close();

    // space steps + output

    fout.open(out_path + "_x1.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
    fout << std::setprecision(8);

    // init space vector
    std::vector<T> x1;

    fout << "x1\n";
    x1.resize(std::size_t( floor((poisson_eq._x1_L1 - poisson_eq._x1_0) / poisson_eq._h1) + 1 ));
    x1[0] = poisson_eq._x1_0;
    fout << x1[0] << "\n";
    for (std::size_t i = 1; i < x1.size() - 1; ++i) {
        x1[i] = x1[i - 1] + poisson_eq._h1;
        fout << x1[i] << "\n";
    }
    x1[x1.size() - 1] = poisson_eq._x1_L1;
    fout << x1[x1.size() - 1] << "\n";
    fout.close();

    std::vector<T> x2;

    fout.open(out_path + "_x2.csv");

    fout << "x2\n";
    x2.resize(std::size_t( floor((poisson_eq._x2_L2 - poisson_eq._x2_0) / poisson_eq._h2) + 1 ));
    x2[0] = poisson_eq._x2_0;
    fout << x2[0] << "\n";
    for (std::size_t i = 1; i < x2.size() - 1; ++i) {
        x2[i] = x2[i - 1] + poisson_eq._h2;
        fout << x2[i] << "\n";
    }
    x2[x2.size() - 1] = poisson_eq._x2_L2;
    fout << x2[x2.size() - 1] << "\n";

    // std::cout << Vector(x) << "\n";

    fout.close();

    // Solution

    fout.open(out_path + "_y.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
    fout << std::setprecision(8);

    // init solution vector
    Matrix<T> y(x1.size(), x2.size());
    std::cout << "\n" << y;
    //output
    for (std::size_t i = 0; i < x1.size(); ++i) {
        for (std::size_t j = 0; j < x2.size(); ++j) {
            if (i == 0 && j == 0){
                fout << "y_" << i << "_" << j;
            }
            else{
                fout << ",y_" << i << "_" << j;
            }
        }
    }
    fout << "\n";

    // get initial values for solution vector + output
    for(std::size_t i = 0; i < x1.size(); ++i){
        y[i][0] = poisson_eq._lower_boundary_condition.second(x1[i], x2[0], time[0]);
        y[i][x2.size() - 1] = poisson_eq._upper_boundary_condition.second(x1[i], x2[x2.size() - 1], time[0]);
    }

    for(std::size_t i = 0; i < x2.size(); ++i){
        y[0][i] = poisson_eq._left_boundary_condition.second(x1[0], x2[i], time[0]);
        y[x1.size() - 1][i] = poisson_eq._right_boundary_condition.second(x1[x1.size() - 1], x2[i], time[0]);
    }
    
    for (std::size_t i = 0; i < x1.size(); ++i) {
        for (std::size_t j = 0; j < x2.size(); ++j) {        
            if(i == 0 && j == 0){
                fout << y[i][j];
            }
            else{
                fout << "," << y[i][j];
            }
        }
    }    
    fout << "\n";

    /**
     * -------------------------------------------------- *
     */

    /**
     *
     *
     * SOLUTION
     *
     *
     */

    poisson_FDM_scheme(poisson_eq, time, x1, x2, y, tolerance, fout);
    
    /**
     * -------------------------------------------------- *
     */

    fout.close();

    return;
};
