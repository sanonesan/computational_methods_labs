#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "../Class_Fred.hpp"
#include "./Fred_schemes.hpp"

template <typename T>
void Fred_solve(
    const Class_Fred<T>& Fred_eq,
    const std::size_t scheme,
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

    if (Fred_eq._Singular == 0){

        // space steps + output

        fout.open(out_path + "_x.csv");
        if (!fout) {
            std::cout << "\n error \n";
            return;
        }

        // init space vector
        std::vector<T> x;

        fout << "x\n";
        x.resize(std::size_t( floor((Fred_eq._x_L - Fred_eq._x_0) / Fred_eq._h) + 1 ));
        x[0] = Fred_eq._x_0;
        fout << x[0] << "\n";
        for (std::size_t i = 1; i < x.size() - 1; ++i) {
            x[i] = x[i - 1] + Fred_eq._h;
            fout << x[i] << "\n";
        }
        x[x.size() - 1] = Fred_eq._x_L;
        fout << x[x.size() - 1] << "\n";
        fout.close();

        std::vector<T> s(x);
        fout.open(out_path + "_s.csv");


        if (!fout) {
            std::cout << "\n error \n";
            return;
        }

        fout << "s\n";

        for (std::size_t i = 0; i < s.size(); ++i) {
            fout << s[i] << "\n";
        }
        fout.close();


        fout.open(out_path + "_y.csv");
        

        if (!fout) {
            std::cout << "\n error \n";
            return;
        }

        fout << "y\n";


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
    
        if (Fred_eq._Invertible == 0){
            if (scheme == 0){
                Fred_non_inv_quadrature_scheme(Fred_eq, x, s, tolerance, fout);
            }
        }

        if (Fred_eq._Invertible == 0){
            if (scheme == 1){
                Fred_non_inv_simple_iter_scheme(Fred_eq, x, s, tolerance, fout);
            }
        }

    } else {


        fout.open(out_path + "_c.csv");
        if (!fout) {
            std::cout << "\n error \n";
            return;
        }
        fout << "c_x,c_y\n";

        for (std::size_t i = 0; i < Fred_eq.c_i.size(); ++i) {
            fout << Fred_eq.c_i[i].first << "," << Fred_eq.c_i[i].second  << "\n";
        }

        fout.close();


        fout.open(out_path + "_y.csv");
        

        if (!fout) {
            std::cout << "\n error \n";
            return;
        }

        fout << "y\n";



        Fred_singular_scheme(Fred_eq, tolerance, fout);
    }



    // poisson_FDM_scheme(Fred_eq, time, x, x2, y, tolerance, fout);
    
    /**
     * -------------------------------------------------- *
     */

    fout.close();

    return;
};
