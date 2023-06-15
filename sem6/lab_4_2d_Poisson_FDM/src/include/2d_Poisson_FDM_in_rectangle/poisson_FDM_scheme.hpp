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

    T norm, old_norm, new_norm, nu;
    


    if (poisson_eq._left_boundary_condition.first == 1) {
        diag_matrix_x1[0][0] = 1.;
        diag_matrix_x1[1][0] = 0;

    }
    else {
        diag_matrix_x1[0][0] = 1. / poisson_eq._tau + 1. / poisson_eq._h1_2;
        diag_matrix_x1[1][0] = - 1. / poisson_eq._h1_2;
    }
    for (std::size_t i = 1; i < x1.size() - 1; ++i) {
        diag_matrix_x1[1][i] = 1. / poisson_eq._h1_2;
        diag_matrix_x1[0][i] = - 2. * ( 1. / poisson_eq._h1_2 + 1. / poisson_eq._tau_2);
        diag_matrix_x1[2][i] = 1. / poisson_eq._h1_2;

    }
    if (poisson_eq._right_boundary_condition.first == 1) {
        diag_matrix_x1[0][diag_matrix_x1[0].size() - 1] = 1.;
        diag_matrix_x1[2][diag_matrix_x1[1].size() - 1] = 0.;

    }
    else {
        diag_matrix_x1[0][diag_matrix_x1[0].size() - 1] = 1. / poisson_eq._tau + 1. / poisson_eq._h1_2;
        diag_matrix_x1[2][diag_matrix_x1[1].size() - 1] = - 1. / poisson_eq._h1_2;
    }




    if (poisson_eq._lower_boundary_condition.first == 1) {
        diag_matrix_x2[0][0] = 1.;
        diag_matrix_x2[1][0] = 0.;
    }
    else {
        diag_matrix_x2[0][0] = 1. / poisson_eq._tau + 1. / poisson_eq._h2_2;
        diag_matrix_x2[1][0] = - 1. / poisson_eq._h2_2;
    }
    for (std::size_t i = 1; i < x2.size() - 1; ++i) {
        diag_matrix_x2[1][i] = 1. / poisson_eq._h2_2;
        diag_matrix_x2[0][i] = - 2. * ( 1. / poisson_eq._h2_2 + 1. / poisson_eq._tau_2);
        diag_matrix_x2[2][i] = 1. / poisson_eq._h2_2;
    }
    if (poisson_eq._upper_boundary_condition.first == 1) {
        diag_matrix_x2[0][diag_matrix_x2[0].size() - 1] = 1.;
        diag_matrix_x2[2][diag_matrix_x2[0].size() - 1] = 0.;
        
    }
    else {
        diag_matrix_x2[0][diag_matrix_x2[0].size() - 1] = 1. / poisson_eq._tau + 1. / poisson_eq._h2_2;
        diag_matrix_x2[2][diag_matrix_x2[0].size() - 1] = - 1. / poisson_eq._h2_2;
    }
    

    std::size_t iter = 0;
    std::size_t k = 0;
    
    while(true){


        for (std::size_t j = 1; j < x2.size() - 1; ++j) {
            
            if (poisson_eq._left_boundary_condition.first == 1)  {
                RHS_1[0] = poisson_eq._left_boundary_condition.second(x1[0], x2[j], time[k]);
            }
            else{
                //что-то если ГУ не 1 рода
                RHS_1[0] = 1. / poisson_eq._tau * y_k[0][j] + 1./ poisson_eq._h1 * poisson_eq._left_boundary_condition.second(x1[0], x2[j], time[k]) 
                                + 1. / 2. / poisson_eq._h2_2 * (y_k[0][j + 1] - 2 * y_k[0][j] + y_k[0][j - 1])  
                                + 1. / 4. * (poisson_eq._f(x1[0], x2[j], time[k]))  + poisson_eq._f(x1[0] + poisson_eq._h1 / 2, x2[j], time[k]);
            }
            
            for (std::size_t i = 1; i < x1.size() - 1; ++i) {
                RHS_1[i] = - (1. / poisson_eq._h2_2 * y_k[i][j + 1] + 2. * (1. / poisson_eq._tau - 1. / poisson_eq._h2_2) * y_k[i][j] 
                            + 1. / poisson_eq._h2_2 * y_k[i][j - 1] + poisson_eq._f(x1[i], x2[j], time[k]));
                
                //- 2 / poisson_eq._tau * y_k[i][j] - (y_k[i][j + 1] - 2 * y_k[i][j] + y_k[i][j - 1]) / poisson_eq._h2_2 - poisson_eq._f(x1[i], x2[j], time[k]);
            }
            
            if (poisson_eq._right_boundary_condition.first == 1)  {
                RHS_1[x1.size() - 1] = poisson_eq._right_boundary_condition.second(x1[x1.size() - 1], x2[j], time[k]);
            }
            else{
                //что-то если ГУ не 1 рода
                RHS_1[x1.size() - 1] = 1. / poisson_eq._tau * y_k[x1.size() - 1][j] + 1. / poisson_eq._h1 * poisson_eq._right_boundary_condition.second(x1[x1.size() - 1], x2[j], time[k]) 
				                + 1. / 2. / poisson_eq._h2_2 * (y_k[x1.size() - 1][j + 1] - 2 * y_k[x1.size() - 1][j] + y_k[x1.size() - 1][j - 1]) 
					            + 1. / 4. * (poisson_eq._f(x1[x1.size() - 1], x2[j], time[k]) + poisson_eq._f(x1[x1.size() - 1] - poisson_eq._h1 / 2, x2[j], time[k]));

			}

            tmp1 = solver.solve_banded(1, 1, diag_matrix_x1, RHS_1); 
            
            for (std::size_t i = 0; i < x1.size(); ++i) {
                y_k_1_2[i][j] = tmp1[i];            
            }   

        }

        for (std::size_t i = 0; i < x1.size(); ++i) {
            y_k_1_2[i][0] = y_k[i][0];
            y_k_1_2[i][x2.size() - 1] = y_k[i][x2.size() - 1];
        }
        


        for (std::size_t i = 1; i < x1.size() - 1; ++i) {
            
            if (poisson_eq._lower_boundary_condition.first == 1)  {
                RHS_2[0] = poisson_eq._lower_boundary_condition.second(x1[i], x2[0], time[k]);
            }
            else{
                //что-то если ГУ не 1 рода
                RHS_2[0] = 1. / poisson_eq._tau * y_k_1_2[i][0] + 1./ poisson_eq._h2 * poisson_eq._lower_boundary_condition.second(x1[i], x2[0], time[k]) 
                                + 1. / 2. / poisson_eq._h1_2 * (y_k_1_2[i + 1][0] - 2 * y_k_1_2[i][0] + y_k_1_2[i - 1][0]) 
                                + 1. / 4. * (poisson_eq._f(x1[i], x2[0], time[k]) + poisson_eq._f(x1[i], x2[0] + poisson_eq._h2 / 2, time[k]));

            }
            
            for (std::size_t j = 1; j < x2.size() - 1; ++j) {
                RHS_2[j] = - (1. / poisson_eq._h1_2 * y_k_1_2[i + 1][j] + 2. * (1. / poisson_eq._tau - 1. / poisson_eq._h1_2) * y_k_1_2[i][j] 
                            + 1. / poisson_eq._h1_2 * y_k_1_2[i - 1][j] + poisson_eq._f(x1[i], x2[j], time[k]));
                ;
                
                //- 2 / poisson_eq._tau * y_k_1_2[i][j] - (y_k_1_2[i + 1][j] - 2 * y_k_1_2[i][j] + y_k_1_2[i - 1][j]) / poisson_eq._h1_2 - poisson_eq._f(x1[i], x2[j], time[k]);
            }
            
            if (poisson_eq._upper_boundary_condition.first == 1)  {
                RHS_2[x2.size() - 1] = poisson_eq._upper_boundary_condition.second(x1[i], x2[x2.size() - 1], time[k]);
            }
            else{
                //что-то если ГУ не 1 рода
                RHS_2[x2.size() - 1] = 1. / poisson_eq._tau * y_k_1_2[i][x2.size() - 1] + 1./ poisson_eq._h2 * poisson_eq._upper_boundary_condition.second(x1[i], x2[x2.size() - 1], time[k]) 
                                + 1. / 2. / poisson_eq._h1_2 * (y_k_1_2[i + 1][x2.size() - 1] - 2 * y_k_1_2[i][x2.size() - 1] + y_k_1_2[i - 1][x2.size() - 1]) 
                                + 1. / 4. * (poisson_eq._f(x1[i], x2[x2.size() - 1], time[k]) + poisson_eq._f(x1[i], x2[x2.size() - 1] - poisson_eq._h2 / 2, time[k]));
            }


            tmp2 = solver.solve_banded(1, 1, diag_matrix_x2, RHS_2); 

            for(std::size_t j = 1; j < x2.size() - 1; ++j){
                y_k_1[i][j] = tmp2[j];
            }

        }


        for (std::size_t j = 0; j < x2.size(); ++j) {
            y_k_1[0][j] = y_k_1_2[0][j];
            y_k_1[x1.size() - 1][j] = y_k_1_2[x1.size() - 1][j];
        }


        norm = 0.;
        
        for (std::size_t i = 0; i < x1.size(); ++i) {
            for (std::size_t j = 0; j < x2.size(); ++j) {
                norm = std::max(norm, fabs(y_k_1[i][j] - y_k[i][j]));
            }
        }




        ++iter;
        ++k;
        std::cout << iter << "\n";
        if (iter != 1){
            old_norm = new_norm;
            new_norm = norm;
            nu = new_norm - old_norm;
        } else {
            new_norm = norm;
        }

        if (iter > 200 || k > time.size() - 2 || new_norm < (1e-16) * (1 - nu)){
            break;
        }


        y_k = y_k_1;

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




};
