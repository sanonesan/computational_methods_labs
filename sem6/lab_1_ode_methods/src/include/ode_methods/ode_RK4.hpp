#ifndef ODE_RK4_HPP
#define ODE_RK4_HPP

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "vector_norm.hpp"

#define eps 1e-16

/*
    template<typename T>
    T RK4_coef(const T t, const T tau, const std::vector<T>& x, const std::size_t k, const F &func);

    template<typename T, typename F>
    void RK4_fix_step(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const string &out_path);

    template<typename T, typename F>
    void RK4_vary_step(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const T tol, const string &out_path);
*/


/*
        Вычисление коэффициентов метода Рунге-Кутты
*/
template <typename T, typename F>
T RK4_coef(const T t, const T tau, const std::vector<T> &x, const std::size_t k, const F &func) {
    std::vector<T> tmp(x);
    T k1, k2, k3, k4;

    k1 = tau * func(tmp, t);

    tmp[k] = x[k] + k1 / 2;
    k2 = tau * func(tmp, t + tau / 2);
    
    tmp[k] = x[k] + k2 / 2;
    k3 = tau * func(tmp, t + tau / 2);
    
    tmp[k] = x[k] + k3;
    k4 = tau * func(tmp, t + tau);

    return (k1 + k4) / 6 + (k2 + k3) / 3;
}

/*
        Метод Рунге-Кутты 4 порядка с фиксированным шагом
*/
template <typename T, typename F>
void ode_RK4_fix_step(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const std::string &out_path) {
    std::ofstream fout(out_path);
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
    fout << std::setprecision(48);
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

    std::vector<T> tmp(x);

    while (start_time <= end_time) {
        fout << start_time + tau;

        for (std::size_t i = 0; i < func.size(); ++i) {
            x[i] += RK4_coef(start_time, tau, tmp, i, func[i]);
            fout << "," << x[i];
        }
        fout << "," << tau << "\n";

        tmp.assign(x.begin(), x.end());

        start_time += tau;
        if(fabs(start_time - end_time) < 1e-8){
            std::cout.precision(17);
            std::cout << start_time << " " << end_time << "\n";

        }
    }

    fout.close();

    return;
}

/*
        Метод Рунге-Кутты 4 порядка с изменяющимя шагом
*/
template <typename T, typename F>
void ode_RK4_vary_step(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const T tol, const std::string &out_path) {
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

    std::vector<T> tmp(x);
    std::vector<T> x_1(x);
    std::vector<T> x_2(x);
    T breaker;
    // T coef_tol = tol / 50000;

    while (true) {
        if (start_time + tau > end_time) {
            if ((start_time < end_time)) {
                tau = end_time - start_time;
                fout << start_time + tau;
                for (std::size_t i = 0; i < func.size(); ++i) {
                    x[i] += RK4_coef(start_time, tau, x_1, i, func[i]);
                    fout << "," << x[i];
                }
                fout << "," << tau << "\n";
            }
            break;
        }

        while (true) {
            for (std::size_t i = 0; i < func.size(); ++i) {
                x_1[i] = x[i] + RK4_coef(start_time, tau, x, i, func[i]);
                tmp[i] = x[i] + RK4_coef(start_time, tau / 2, x, i, func[i]);
            }

            for (std::size_t i = 0; i < func.size(); ++i)
                x_2[i] = tmp[i] + RK4_coef(start_time + tau / 2, tau / 2, tmp, i, func[i]);

            for (std::size_t i = 0; i < func.size(); ++i) {
                tmp[i] = x_1[i] - x_2[i];
            }

            breaker = norm_inf(tmp) / 15;
            if (breaker >= tol) {
                tau /= 2;
            } else {
                x.assign(x_1.begin(), x_1.end());
                start_time += tau;

                fout << start_time;
                for (std::size_t i = 0; i < func.size(); ++i) {
                    fout << "," << x[i];
                }
                fout << "," << tau << "\n";

                tau *= 2;

                break;
            }
        }

        
    }

    fout.close();

    return;
}

#endif