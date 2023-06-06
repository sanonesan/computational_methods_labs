#pragma once

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "../Class_2d_Poisson_equation_in_rectangle.hpp"

template <typename T>
void poisson_FDM_scheme(

    const Class_2d_Poisson_equation_in_rectangle<T>& poisson_eq,
    const std::vector<T>& time,
    const std::vector<T>& x1,
    const std::vector<T>& x2,
    Matrix<T> y_k,
    const T& tolerance,
    std::ofstream& fout) {
    Solver_SLE<T> solver;

    //solver.solve_banded();

    Matrix<T> diag_matrix_x1(3, x1.size());
    Matrix<T> diag_matrix_x2(3, x2.size());

    Vector<T> RHS_1(x1.size());
    Vector<T> RHS_2(x2.size());
    Vector<T> tmp1(x1.size());
    Vector<T> tmp2(x2.size());
    

    Matrix<T> y_k_1_2(y_k);
    Matrix<T> y_k_1(y_k);

    diag_matrix_x1[0][0] = 1.;
    for (std::size_t i = 1; i < x1.size() - 1; ++i) {
        diag_matrix_x1[1][i] = 1. / poisson_eq._h1_2;
        diag_matrix_x1[0][i] = - 2. * ( 1. / poisson_eq._h1_2 + 1. / poisson_eq._tau_2);
        diag_matrix_x1[2][i] = 1. / poisson_eq._h1_2;

    }
    diag_matrix_x1[0][diag_matrix_x1[0].size() - 1] = 1.;

    diag_matrix_x2[0][0] = 1.;
    for (std::size_t i = 1; i < x2.size() - 1; ++i) {
        diag_matrix_x2[1][i] = 1. / poisson_eq._h2_2;
        diag_matrix_x2[0][i] = - 2. * ( 1. / poisson_eq._h2_2 + 1. / poisson_eq._tau_2);
        diag_matrix_x2[2][i] = 1. / poisson_eq._h2_2;
    }
    diag_matrix_x2[0][diag_matrix_x2[0].size() - 1] = 1.;


    for (std::size_t k = 1; k < time.size(); ++k){
        
        //нижнее ГУ 1 рода
        if (poisson_eq._lower_boundary_condition.first == 1){
            for (std::size_t i = 0; i < x1.size(); ++i) {
                y_k_1_2[i][0] = poisson_eq._lower_boundary_condition.second(x1[i], x2[0], time[k]);                
            }
        }
        else {
            //что-то если ГУ не 1 рода
        }

        for (std::size_t j = 1; j < x2.size() - 1; ++j) {
            
            if (poisson_eq._left_boundary_condition.first == 1)  {
                RHS_1[0] = poisson_eq._left_boundary_condition.second(x1[0], x2[j], time[k]);
            }
            else{
                //что-то если ГУ не 1 рода
            }
            
            for (std::size_t i = 1; i < x1.size() - 1; ++i) {
                RHS_1[i] = - 2 / poisson_eq._tau * y_k[i][j] - (y_k[i][j + 1] - 2 * y_k[i][j] + y_k[i][j - 1]) / poisson_eq._h2_2 - poisson_eq._f(x1[i], x2[j], time[k]);
            }
            
            if (poisson_eq._right_boundary_condition.first == 1)  {
                RHS_1[x1.size() - 1] = poisson_eq._right_boundary_condition.second(x1[x1.size() - 1], x2[j], time[k]);
            }
            else{
                //что-то если ГУ не 1 рода
            }

            tmp1 = solver.solve_banded(1, 1, diag_matrix_x1, RHS_1);    
            for (std::size_t i = 0; i < x1.size(); ++i) {
                y_k_1_2[i][j] = tmp1[i];            
            }   


        }

        //верхнее ГУ 1 рода
        if (poisson_eq._upper_boundary_condition.first == 1){
            for (std::size_t i = 0; i < x1.size(); ++i) {
                y_k_1_2[i][x2.size() - 1] = poisson_eq._upper_boundary_condition.second(x1[i], x2[x2.size() - 1], time[k]);                
            }
        }
        else {
            //что-то если ГУ не 1 рода
        }


        //левое ГУ 1 рода
        if (poisson_eq._left_boundary_condition.first == 1){
            for (std::size_t j = 0; j < x2.size(); ++j) {
                y_k_1[0][j] = poisson_eq._left_boundary_condition.second(x1[0], x2[j], time[k]);                
            }
        }
        else {
            //что-то если ГУ не 1 рода
        }

        for (std::size_t i = 1; i < x1.size() - 1; ++i) {
            
            if (poisson_eq._lower_boundary_condition.first == 1)  {
                RHS_2[0] = poisson_eq._lower_boundary_condition.second(x1[i], x2[0], time[k]);
            }
            else{
                //что-то если ГУ не 1 рода
            }
            
            for (std::size_t j = 1; j < x2.size() - 1; ++j) {
                RHS_2[j] = - 2 / poisson_eq._tau * y_k_1_2[i][j] - (y_k_1_2[i + 1][j] - 2 * y_k_1_2[i][j] + y_k_1_2[i - 1][j]) / poisson_eq._h1_2 - poisson_eq._f(x1[i], x2[j], time[k]);
            }
            
            if (poisson_eq._upper_boundary_condition.first == 1)  {
                RHS_2[x2.size() - 1] = poisson_eq._upper_boundary_condition.second(x1[i], x2[x2.size() - 1], time[k]);
            }
            else{
                //что-то если ГУ не 1 рода
            }


            tmp2 = solver.solve_banded(1, 1, diag_matrix_x2, RHS_2); 

            for(std::size_t j = 1; j < x2.size() - 1; ++j){
                y_k_1[i][j] = tmp2[j];
            }


        }


        //правое ГУ 1 рода
        if (poisson_eq._right_boundary_condition.first == 1){
            for (std::size_t j = 0; j < x2.size(); ++j) {
                y_k_1[x1.size() - 1][j] = poisson_eq._right_boundary_condition.second(x1[x1.size() - 1], x2[j], time[k]);                
            }
        }
        else {
            //что-то если ГУ не 1 рода
        }


        for (std::size_t i = 0; i < x1.size(); ++i) {
            for (std::size_t j = 0; j < x2.size(); ++j) {        
                if(i == 0 && j == 0){
                    fout << y_k_1[i][j];
                }
                else{
                    fout << "," << y_k_1[i][j];
                }
            }
        }    
        fout << "\n";

        y_k = y_k_1;

    }
        


        
};
