#ifndef SOLVER_SLE
#define SOLVER_SLE

#include <tuple>
#include <vector>

#include "../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "./SLE_methods/QR_method.hpp"
#include "./SLE_methods/method_gauss.hpp"

/**
 * Решение СЛАУ:
 * - метод Гаусса;
 * - метод QR-разложения.
 **/
template <class T>
class Solver_SLE {

public:
    double tol = 1e-3;
    std::string file_name = "";
    std::string output_folder = "../output/";

    // Default constructor
    Solver_SLE() {
        this->tol = 1e-9;
        this->output_folder = "../output/";
        this->file_name = "";
    }

    // Alt constructor
    Solver_SLE(T tol, std::string output_folder, std::string file_name) {
        this->tol = tol;
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

};

#endif