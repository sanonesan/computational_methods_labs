#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "../Class_1d_heat_equation.hpp"
#include "./explicit_2_layer_difference_scheme.hpp"
#include "./implicit_2_layer_difference_scheme.hpp"
#include "./mixed_2_layer_difference_scheme.hpp"

template <typename T, typename sigma_type>
void templated_2_layer_difference_scheme(
    const Class_1d_heat_equation<T>& heat_equation,
    const sigma_type& sigma,
    const std::string& out_path) {
    /**
     *
     *
     * GETTING INITIAL VALUES
     *
     *
     */

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
    time.resize(std::size_t((heat_equation._end_time - heat_equation._start_time) / heat_equation._tau + 1));
    time[0] = heat_equation._start_time;
    fout << time[0] << "\n";
    for (std::size_t i = 1; i < time.size() - 1; ++i) {
        time[i] = time[i - 1] + heat_equation._tau;
        fout << time[i] << "\n";
    }
    time[time.size() - 1] = heat_equation._end_time;
    fout << time[time.size() - 1] << "\n";

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
    x.resize(std::size_t((heat_equation._xL - heat_equation._x0) / heat_equation._h + 1));
    x[0] = heat_equation._x0;
    fout << x[0] << "\n";
    for (std::size_t i = 1; i < x.size() - 1; ++i) {
        x[i] = x[i - 1] + heat_equation._h;
        fout << x[i] << "\n";
    }
    x[x.size() - 1] = heat_equation._xL;
    fout << x[x.size() - 1] << "\n";

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

    if (fabs(sigma) < 1e-16) {
        // std::cout << "explicit scheme \n";
        explicit_2_layer_difference_scheme(heat_equation, time, x, y, fout);

    } else if (fabs(sigma - 1.) < 1e-16) {
        // std::cout << "implicit scheme \n";
        implicit_2_layer_difference_scheme(heat_equation, time, x, y, fout);
    } else {
        // std::cout << "mixed scheme \n";
        mixed_2_layer_difference_scheme(heat_equation, time, x, y, sigma, fout);
    }

    /**
     * -------------------------------------------------- *
     */

    fout.close();

    return;
};
