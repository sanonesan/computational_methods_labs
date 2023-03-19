#ifndef ODE_SYMMETRICAL_SCHEME_HPP
#define ODE_SYMMETRICAL_SCHEME_HPP

#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../../../sem5/lab_1_SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "./numerical_Jacobian_for_ode_sys.hpp"

#define eps 1e-16



/*
    Двухшаговая симметричная схема 2 порядка (предиктор-корректор)
*/
template <typename T, typename F>
void ode_2_step_symmetrical_scheme(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const std::string &out_path, const T &tol) {
    std::ofstream fout(out_path);
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
	fout << std::setprecision(8);
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
    Vector<T> tmp1(x);

    Vector<T> tmp_foo(func.size());

    while (start_time <= end_time) {

        

        //predictor (explicit Euler)
        for (std::size_t i = 0; i < x.size(); ++i) {
            tmp_foo[i] = func[i](x, start_time);
            tmp[i] = x[i] + tau * tmp_foo[i];
        }

        start_time += tau;

        fout << start_time;

        //corrector
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] = x[i] + tau / 2 * (tmp_foo[i] + func[i](tmp, start_time));
            fout << "," << x[i];
        }
        fout << "," << tau << "\n";

    }

    return;
}

/*
    Симметричная схема 2 порядка (решение нелинейного уравнения)
*/
template <typename T, typename F>
void ode_symmetrical_scheme_nonlinear_eq(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const std::string &out_path, const T &tol) {
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
        Jacobian *= -(tau / 2);
        // E - tau/2 * J
        for (std::size_t i = 0; i < Jacobian.get_rows(); ++i) {
            Jacobian[i][i] = 1 + Jacobian[i][i];
        }

        for (std::size_t i = 0; i < func.size(); ++i) {
            tmp[i] = x[i] + tau / 2 * func[i](x, start_time);
        }

        /**
         * Solving equation:
         * (E - tau / 2 * J) * y_{n+1} = F_{n}
         * for y_{n+1}, where
         * F_{n} = y_{n} + (tau/2) * f(t_{n}, y_{n})
         **/
        tmp = std::get<0>(solver.QR(Jacobian, tmp));
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