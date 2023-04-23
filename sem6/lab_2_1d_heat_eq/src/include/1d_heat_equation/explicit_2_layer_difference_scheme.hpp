#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../Class_1d_heat_equation.hpp"
#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"


//sigma == 0
template <typename T>
void explicit_2_layer_difference_scheme(

    const Class_1d_heat_equation<T> &heat_equation,
    const std::vector<T> &x,
    const std::vector<T> &time,
    std::vector<T> &y,
    std::ofstream &fout){

    auto _K_approximation = [heat_equation, y, x](const std::size_t i) -> T{
        if(i == 0){
            throw std::invalid_argument("0 < i < y.size()");
        }
        return (heat_equation._K(y[i], x[i]) + heat_equation._K(y[i-1], x[i-1])) / 2;
    };

    Solver_SLE<T> solver_SLE;
    Vector<T> solution(y.size());
    Vector<T> b(y.size());
    Matrix<T> banded_matrix(1, y.size());

    T c_rho_h_tau = 0.;

    c_rho_h_tau = heat_equation._c * heat_equation._rho * heat_equation._h / heat_equation._tau;

    T a_i = 0.;
    T a_i_1 = 0.;
    T kappa = 0.;
    T mu = 0.;

    for (std::size_t j = 1; j < time.size(); ++j) {
        b.assign(y.begin(), y.end());

        banded_matrix[0][0] = 1;

        if (heat_equation._left_boundary_condition_type == 1) {

            mu = c_rho_h_tau / 2 * y[0];

            a_i = _K_approximation(1);

            mu -= (heat_equation._boundary_conditions[0](x[0], time[j-1]) - a_i * (y[1] - y[0]) / heat_equation._h );
            
            mu /= c_rho_h_tau / 2;
            b[0] = mu;        
            
        }

        for (std::size_t i = 1; i < b.size() - 1; ++i) {

            a_i = _K_approximation(i);
            a_i_1 = _K_approximation(i + 1);

            // diag (-C_i)
            banded_matrix[0][i] = -c_rho_h_tau;

            b[i] = -(c_rho_h_tau * y[i] + (a_i_1 * (y[i + 1] - y[i]) - a_i * (y[i] - y[i - 1])) / heat_equation._h);
        }

        banded_matrix[0][banded_matrix[0].size() - 1] = 1;

        if (heat_equation._right_boundary_condition_type == 1) {

            mu = c_rho_h_tau / 2 * y[y.size() - 1];
            
            a_i = _K_approximation(y.size() - 1);

            mu += (heat_equation._boundary_conditions[1](x[x.size() - 1], time[j-1]) - a_i * (y[y.size() - 1] - y[y.size() - 2]) / heat_equation._h);
            mu /= c_rho_h_tau / 2;
            b[b.size() - 1] = mu;        
            
        }

        solution = solver_SLE.solve_banded(0, 0, banded_matrix, b);

        for (std::size_t i = 0; i < solution.size(); ++i) {
            y[i] = solution[i];
        }

        fout << y[0];
        for (std::size_t i = 1; i < x.size() - 1; ++i) {
            fout << "," << y[i];
        }
        fout << "," << y[y.size() - 1] << "\n";
    }


};
