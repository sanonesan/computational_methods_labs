#ifndef ODE_EXPLICIT_EULER_HPP
#define ODE_EXPLICIT_EULER_HPP

#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <vector>

#define eps 1e-16

/*
    Явный метод Эйлера
*/
template <typename T, typename F>
void ode_explicit_Euler(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const std::string &out_path) {
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

    while (start_time <= end_time) {
        fout << start_time + tau;
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] = tmp[i] + tau * func[i](tmp, start_time);
            fout << "," << x[i];
        }
        fout << "," << tau << "\n";
        tmp.assign(x.begin(), x.end());
        start_time += tau;
    }

    return;
}

#endif
