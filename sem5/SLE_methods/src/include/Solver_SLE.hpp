#pragma once

#ifndef SOLVER_SLE
#define SOLVER_SLE

#include <tuple>
#include <vector>

#include "../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "./SLE_methods/QR_method.hpp"
#include "./SLE_methods/method_gauss.hpp"
#include "./SLE_methods/banded_sle.hpp"


#include "./SLE_methods/Jacoby_method.hpp"
#include "./SLE_methods/simple_iter_method.hpp"

/**
 * Решение СЛАУ:
 * - метод Гаусса;
 * - метод QR-разложения.
 **/
template <class T>
class Solver_SLE {

public:
    T tolerance = 1e-9;
    std::string file_name = "";
    std::string output_folder = "../output/";

    // Default constructor
    Solver_SLE() {
        this->tolerance = 1e-9;
        this->output_folder = "../output/";
        this->file_name = "";
    }

    // Alt constructor
    Solver_SLE(T tolerance, std::string output_folder, std::string file_name) {
        this->tolerance = tolerance;
        this->output_folder = "../output/";
        this->file_name = file_name;
    }

    /**
     * @brief Решение СЛАУ методом Гаусса
     *
     * @result Vector<T> solution
     **/
    Vector<T> Gauss(const Matrix<T> &A, const Vector<T> &b) {
        // std::string out_path;
        // out_path = this->output_folder + this->file_name + ".csv";

        Vector<T> result;

        method_gauss(A, b, result);

        return result;
    }

    /**
     * @brief Решение СЛАУ методом QR-разложения
     *
     * @result std::tuple<Vector<T> solution, Matrix<T> Q, Matrix<T> R>
     **/
    std::tuple<Vector<T>, Matrix<T>, Matrix<T>> QR(const Matrix<T> &A, const Vector<T> &b) {
        // std::string out_path;
        // out_path = this->output_folder + this->file_name + ".csv";

        Matrix<T> Q(A.get_rows(), A.get_cols());
        Matrix<T> R(A.get_rows(), A.get_cols());
        Vector<T> solution;

        QR_method(A, b, Q, R, solution);

        std::tuple<Vector<T>, Matrix<T>, Matrix<T>> result(solution, Q, R);

        return result;
    }

    /**
     * @brief Get Q & R decomposion
     *
     * @result std::tuple<solution, Matrix<T> Q, Matrix<T> R>
     **/
    std::tuple<Matrix<T>, Matrix<T>> QR_decomposion(const Matrix<T> &A) {
        // std::string out_path;
        // out_path = this->output_folder + this->file_name + ".csv";

        Matrix<T> Q(A.get_rows(), A.get_cols());
        Matrix<T> R(A.get_rows(), A.get_cols());
        Vector<T> solution;

        QR_decomposion_method(A, Q, R);

        std::tuple<Matrix<T>, Matrix<T>> result(Q, R);

        return result;
    }

    /**   Ax = b
    *
    *    Regular matrix A: 
    *    [b, c, 0,........]   [d]
    *    [a, b, c,........]   [d]
    *    [0, a, b, c,.....]   [d]
    *    .................. = ...
    *    [......., a, b, c]   [d]
    *    [.........., a, b]   [d]  
    * 
    *    Vector b:
    * 
    *    b = [d, d, d, d, d]
    * 
    * 
    *   Banded matrix from A:
    * 
    *    banded_matrix = [
    *        [c, c, c, c, *],
    *        [b, b, b, b, b],
    *        [*, a, a, a, a],
    *    ]
    *    
    */
    Vector<T> solve_banded(const std::size_t l, const std::size_t u, Matrix<T>& banded_matrix, Vector<T>& b){
        
        Vector<T> solution = banded_sle(l, u, banded_matrix, b);

        return solution;
    };


    Vector<T> solve_Jacoby(Matrix<T>& A, Vector<T>& b, Vector<T>& x0){
        
        Vector<T> solution = Jacoby_method(A, b, x0, this->tolerance);

        return solution;
    };

    Vector<T> solve_simple_iter(Matrix<T>& A, Vector<T>& b, Vector<T>& x0, const T tau, const std::size_t M){
        
        Vector<T> solution = simple_iter_method(A, b, x0, tau, M, this->tolerance);

        return solution;
    };

};

#endif