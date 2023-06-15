#pragma once

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "../Class_Fred.hpp"

template <typename T>
void Fred_non_inv_quadrature_scheme(

    const Class_Fred<T>& Fred_eq,
    const std::vector<T>& x,
    const std::vector<T>& s,
    const T& tolerance,
    std::ofstream& fout) {


    Solver_SLE<T> solver_SLE;


    //Solution
    Vector<T> y(x);

    Matrix<T> Sys(x.size(), x.size());
    Sys.make_matrix_identity(x.size());
    Vector<T> RHS(x.size());

    Vector<T> coefs(x);

    /**
     * Формула трапеций
    */
    coefs[0] = Fred_eq._h;
    for (std::size_t i = 1; i < coefs.size() - 1; ++i) {
        coefs[i] = Fred_eq._h / 2;
    }
    coefs[coefs.size() - 1] = Fred_eq._h;



    for (std::size_t i = 0; i < x.size(); ++i) {

        Sys[i][0] -= coefs[0] * Fred_eq._lambda * Fred_eq._K(x[i], s[0]);
        for (std::size_t j = 1; j < x.size()-1; ++j) {
            Sys[i][j] -= coefs[j] * Fred_eq._lambda * (Fred_eq._K(x[i], s[j]) + Fred_eq._K(x[i], s[j - 1]));
        }
        Sys[i][coefs.size() - 1] -= coefs[coefs.size() - 1] * Fred_eq._lambda * Fred_eq._K(x[i], s[coefs.size() - 1]);

        RHS[i] = Fred_eq._f(x[i]);
    }

    y = solver_SLE.Gauss(Sys, RHS);

    for (std::size_t i = 0; i < y.size(); ++i) {
        fout << y[i] << "\n";
    }


};



template <typename T>
void Fred_non_inv_simple_iter_scheme(

    const Class_Fred<T>& Fred_eq,
    const std::vector<T>& x,
    const std::vector<T>& s,
    const T& tolerance,
    std::ofstream& fout) {


    Solver_SLE<T> solver_SLE;


    //Solution
    Vector<T> y(x);
    Vector<T> y_next(x);

    Vector<T> coefs(x);

    
    // /**
    //  * Правые прямоугольники
    // */
    // for (std::size_t i = 0; i < coefs.size(); ++i) {
    //     coefs[i] = Fred_eq._h;
    //     y[i] = 1;
    // }


    /**
     * Формула трапеций
    */
    coefs[0] = Fred_eq._h;
    for (std::size_t i = 1; i < coefs.size() - 1; ++i) {
        coefs[i] = Fred_eq._h / 2;
    }
    coefs[coefs.size() - 1] = Fred_eq._h;




    while (true)
    {

        for (std::size_t i = 0; i < x.size(); ++i) {
            y_next[i] = Fred_eq._f(x[i]);

            y_next[i] += coefs[0] * Fred_eq._lambda * Fred_eq._K(x[i], s[0]) * y[0];

            for (std::size_t j = 1; j < x.size() - 1; ++j) {
                y_next[i] += coefs[j] * Fred_eq._lambda * (Fred_eq._K(x[i], s[j]) + Fred_eq._K(x[i], s[j])) * y[j];
            }

            y_next[i] += coefs[coefs.size() - 1] * Fred_eq._lambda * Fred_eq._K(x[i], s[coefs.size() - 1]) * y[coefs.size() - 1];
        }


        if ( (y - y_next).norm_inf() < 1e-9){
            y = y_next;
            break;
        }
        y = y_next;
    }


    for (std::size_t i = 0; i < y.size(); ++i) {
        fout << y[i] << "\n";
    }


};


template <typename T>
void Fred_singular_scheme(

    const Class_Fred<T>& Fred_eq,
    const T& tolerance,
    std::ofstream& fout) {


    Solver_SLE<T> solver_SLE;

    //Solution
    Vector<T> y(Fred_eq.c_i.size() + 1);

    Matrix<T> Sys(Fred_eq.c_i.size() + 1, Fred_eq.c_i.size() + 1);
    Sys.make_matrix_identity(Fred_eq.c_i.size() + 1);
    Vector<T> RHS(Fred_eq.c_i.size() + 1);

    std::pair<T,T> tmp;

    for (std::size_t i = 0; i < Fred_eq.k_i.size(); ++i) {


        for (std::size_t j = 0; j < Fred_eq.c_i.size(); ++j) {
            tmp = Fred_eq._Q(Fred_eq.k_i[i], Fred_eq.c_i[j]);
            Sys[i][j] = (Fred_eq.n_i[i].first * tmp.first 
                        + Fred_eq.n_i[i].second  * tmp.second) * Fred_eq.delta_l;    
        }

        Sys[i][Fred_eq.c_i.size()] = 1.;

        RHS[i] = Fred_eq._f_singular(Fred_eq.k_i[i]);

    }
    for (std::size_t j = 0; j < Fred_eq.c_i.size(); ++j) {
        Sys[Fred_eq.c_i.size()][j] = Fred_eq.delta_l;
    }

    y = solver_SLE.Gauss(Sys, RHS);

    for (std::size_t i = 0; i < y.size() - 1; ++i) {
        fout << y[i] << "\n";
    }

    // std::cout << std::scientific;
    // std::cout << std::setprecision(32);
    // std::cout << "R = " << y[y.size() - 1] << "\n";
    


};
