#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../Class_1d_heat_equation.hpp"
#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"


//0 < sigma < 1
template <typename T, typename sigma_type>
void mixed_2_layer_difference_scheme(

    const Class_1d_heat_equation<T> &heat_equation,
    const std::vector<T> &time,
    const std::vector<T> &x,
    std::vector<T> &y,
    const sigma_type &sigma,
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
    Matrix<T> banded_matrix(3, y.size());

    T sigma_h = 1.;
    T c_rho_h_tau = 0.;

    sigma_h /= heat_equation._h;
    c_rho_h_tau = heat_equation._c * heat_equation._rho * heat_equation._h / heat_equation._tau;

    T a_i = 0.;
    T a_i_1 = 0.;
    T kappa = 0.;
    T mu = 0.;

    for (std::size_t j = 1; j < time.size(); ++j) {
        
        banded_matrix[0][0] = 1;

        if (heat_equation._left_boundary_condition_type == 0) {
            b[0] = heat_equation._boundary_conditions[0](x[0], time[j - 1]);
        } else if (heat_equation._left_boundary_condition_type == 1) {
            a_i = _K_approximation(1);

            mu = c_rho_h_tau / 2 * y[0];
            mu -= sigma * heat_equation._boundary_conditions[0](x[0], time[j]);
            mu -= (1 - sigma) * (heat_equation._boundary_conditions[0](x[0], time[j-1]) - a_i * (y[1] - y[0]) / heat_equation._h );
            mu /= c_rho_h_tau / 2 + sigma_h * a_i;
            b[0] = mu;        
            
            kappa = sigma_h * a_i;
            kappa /= c_rho_h_tau / 2 + sigma_h * a_i;
            banded_matrix[1][0] = -kappa;

        }

        for (std::size_t i = 1; i < b.size() - 1; ++i) {

            a_i = _K_approximation(i);
            a_i_1 = _K_approximation(i + 1);

            // upper (B_i)
            banded_matrix[1][i] = sigma_h * a_i_1;

            // diag (-C_i)
            banded_matrix[0][i] = -(sigma_h * (a_i + a_i_1) + c_rho_h_tau);

            // lower (A_i)
            banded_matrix[2][i] = sigma_h * a_i;

            b[i] = -(c_rho_h_tau * y[i] + (1 - sigma) * (a_i_1 * (y[i + 1] - y[i]) - a_i * (y[i] - y[i - 1])) / heat_equation._h);
        }

        banded_matrix[0][banded_matrix[0].size() - 1] = 1;

        if (heat_equation._right_boundary_condition_type == 0) {
            b[b.size() - 1] = heat_equation._boundary_conditions[1](x[x.size() - 1], time[j - 1]);
        } else if (heat_equation._right_boundary_condition_type == 1) {
            a_i = _K_approximation(y.size() - 1);


            mu = c_rho_h_tau / 2 * y[y.size() - 1];
            mu += sigma * heat_equation._boundary_conditions[1](x[x.size() - 1], time[j]);
            mu += (1 - sigma) * (heat_equation._boundary_conditions[1](x[x.size() - 1], time[j-1]) - a_i * (y[y.size() - 1] - y[y.size() - 2]) / heat_equation._h);
            mu /= c_rho_h_tau / 2 + sigma_h * a_i;
            b[b.size() - 1] = mu;        
            
            kappa = sigma_h * a_i;
            kappa /= c_rho_h_tau / 2 + sigma_h * a_i;
            banded_matrix[2][banded_matrix[2].size() - 1] = -kappa;

        }

        solution = solver_SLE.solve_banded(1, 1, banded_matrix, b);

        y.assign(solution.begin(), solution.end());

        fout << y[0];
        for (std::size_t i = 1; i < x.size() - 1; ++i) {
            fout << "," << y[i];
        }
        fout << "," << y[y.size() - 1] << "\n";
    }

};