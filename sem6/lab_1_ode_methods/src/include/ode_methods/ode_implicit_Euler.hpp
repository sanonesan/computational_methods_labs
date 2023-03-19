#ifndef ODE_IMPLICIT_EULER_HPP
#define ODE_IMPLICIT_EULER_HPP

#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

#include "../../../../../sem5/lab_1_SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "./numerical_Jacobian_for_ode_sys.hpp"

#define eps 1e-16


/*
    Неявный метод Эйлера
*/
template <typename T, typename F>
void ode_implicit_Euler(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const std::string &out_path, const T &tol) {
    std::ofstream fout(out_path);
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
    fout << "time";
    for (std::size_t i = 0; i < func.size(); ++i) {
        fout << ",u" << i;
    }
    fout << ",tau\n";
    fout << start_time;
    for (std::size_t i = 0; i < func.size(); ++i) {
        fout << "," << x[i];
    }
    fout << "," << tau << "\n";

    Vector<T> tmp(x);
    Solver_SLE<T> solver;
    Matrix<T> Jacobian;

    while (start_time <= end_time) {
        fout << start_time + tau;

        // J
        Jacobian = numerical_Jacobian_for_ode_sys(func, x, start_time);
        // tau * J
        Jacobian *= -tau;
        // E - tau * J
        for (std::size_t i = 0; i < Jacobian.get_rows(); ++i)
            Jacobian[i][i] = 1 + Jacobian[i][i];

        /**
         * Solving equation:
         * (E - tau * J) * y_{n+1} = y_{n}
         * for y_{n+1}
         **/
        tmp = std::get<0>(solver.QR(Jacobian, x)); 
        x.assign(tmp.begin(), tmp.end());

        for (std::size_t i = 0; i < x.size(); ++i) {
            fout << "," << x[i];
        }
        fout << "," << tau << "\n";

        start_time += tau;
    }

    return;
}

#endif
