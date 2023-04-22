#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "../Class_1d_heat_equation.hpp"

template <typename T>
void templated_2_layer_difference_scheme(const Class_1d_heat_equation<T>& heat_equation, const T sigma, const std::string& out_path) {
    // time steps + output

    std::ofstream fout(out_path + "_t.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
    fout << std::setprecision(8);

    // init time vector
    std::vector<T> time;

    fout << "time\n";
    time.push_back(heat_equation._start_time);
    fout << time[time.size() - 1] << "\n";
    while (time[time.size() - 1] <= heat_equation._end_time) {
        time.push_back(time[time.size() - 1] + heat_equation._tau);
        fout << time[time.size() - 1] << "\n";
    }

    fout.close();

    // space steps + output

    fout.open(out_path + "_x.csv");
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
    fout << std::setprecision(8);

    // init space vector
    std::vector<T> x;

    fout << "x\n";
    x.push_back(heat_equation._x0);
    fout << x[x.size() - 1] << "\n";
    while (x[x.size() - 1] <= heat_equation._xL) {
        x.push_back(x[x.size() - 1] + heat_equation._h);
        fout << x[x.size() - 1] << "\n";
    }

    //std::cout << Vector(x) << "\n";

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
    std::vector<T> y;
    y.assign(x.begin(), x.end());

    // output
    fout << "t0";
    for (std::size_t i = 1; i < x.size(); ++i) {
        fout << ",t" << i;
    }
    fout << "\n";

    // get initial values for solution vector + output
    y[0] = heat_equation._initial_conditions(x[0], time[0]);
    fout << y[0];
    for (std::size_t i = 1; i < x.size() - 1; ++i) {
        y[i] = heat_equation._initial_conditions(x[i], time[0]);
        fout << "," << y[i];
    }
    y[y.size() - 1] = heat_equation._initial_conditions(x[x.size() - 1], time[0]);
    fout << "," << y[y.size() - 1] << "\n";

    // vector for solutions from previous time step
    std::vector<T> y_prev;
    y_prev.assign(y.begin(), y.end());

    // //vector for K(u, x)
    // std::vector<T> a_i;
    // a_i.assign(x.begin(), x.end());
    // for(std::size_t i = 1; i < x.size(); ++i){
    //     a_i[i] = (heat_equation._K(x[i]) + heat_equation._K(x[i-1])) / 2;
    // }

    T sigma_h = 0.;
    T c_rho_h_tau = 0.;

    sigma_h = sigma / heat_equation._h;
    c_rho_h_tau = heat_equation._c * heat_equation._rho * heat_equation._h / heat_equation._tau;

    Solver_SLE<T> solver_SLE;
    Vector<T> solution(y.size());
    Vector<T> b(y.size());
    Matrix<T> banded_matrix(3, y.size());

    // computating solutions + output
    T a_i = 0.;
    T a_i_1 = 0.;
    T kappa = 0.;
    T mu = 0.;

    for (std::size_t j = 1; j < time.size(); ++j) {
        b.assign(y.begin(), y.end());

        banded_matrix[1][0] = 1;

        if (heat_equation._left_boundary_condition_type == 1) {

            mu = c_rho_h_tau / 2 * y[0];
            mu -= sigma * heat_equation._boundary_conditions[0](x[0], time[j]);

            a_i = (heat_equation._K(y[1], x[1]) + heat_equation._K(y[0], x[0])) / 2;

            mu -= (1 - sigma) * (heat_equation._boundary_conditions[0](x[0], time[j-1]) - a_i * (y[1] - y[0]) / heat_equation._h );
            
            mu /= c_rho_h_tau / 2 + sigma_h * a_i;
            b[0] = mu;        
            
            kappa = sigma_h * a_i;
            kappa /= c_rho_h_tau / 2 + sigma_h * a_i;
            banded_matrix[0][0] = -kappa;


        }

        for (std::size_t i = 1; i < b.size() - 1; ++i) {
            a_i = (heat_equation._K(y[i], x[i]) + heat_equation._K(y[i - 1], x[i - 1])) / 2;
            a_i_1 = (heat_equation._K(y[i + 1], x[i + 1]) + heat_equation._K(y[i], x[i])) / 2;

            // upper (B_i)
            if (i < y.size() - 1)
                banded_matrix[0][i] = sigma_h * a_i_1;

            // diag (-C_i)
            banded_matrix[1][i] = -(sigma_h * (a_i + a_i_1) + c_rho_h_tau);

            // lower (A_i)
            if (i > 0)
                banded_matrix[2][i] = sigma_h * a_i;

            b[i] = -(c_rho_h_tau * y[i] + (1 - sigma) * (a_i_1 * (y[i + 1] - y[i]) - a_i * (y[i] - y[i - 1])) / heat_equation._h);
        }

        banded_matrix[1][banded_matrix[1].size() - 1] = 1;

        if (heat_equation._right_boundary_condition_type == 1) {

            mu = c_rho_h_tau / 2 * y[y.size() - 1];
            mu += sigma * heat_equation._boundary_conditions[1](x[x.size() - 1], time[j]);

            a_i = (heat_equation._K(y[y.size() - 1], x[y.size() - 1]) + heat_equation._K(y[y.size() - 2], x[y.size() - 2])) / 2;

            mu += (1 - sigma) * (heat_equation._boundary_conditions[1](x[x.size() - 1], time[j-1]) - a_i * (y[y.size() - 1] - y[y.size() - 2]) / heat_equation._h);
            mu /= c_rho_h_tau / 2 + sigma_h * a_i;
            b[b.size() - 1] = mu;        
            
            kappa = sigma_h * a_i;
            kappa /= c_rho_h_tau / 2 + sigma_h * a_i;
            banded_matrix[2][banded_matrix[2].size() - 1] = -kappa;

        }

        solution = solver_SLE.solve_banded(1, 1, banded_matrix, b);

        for (std::size_t i = 0; i < solution.size(); ++i) {
            y[i] = solution[i];
        }

        fout << y[0];
        for (std::size_t i = 1; i < x.size() - 1; ++i) {
            fout << "," << y[i];
        }
        fout << "," << y[y.size() - 1] << "\n";
    }

    fout.close();
};
