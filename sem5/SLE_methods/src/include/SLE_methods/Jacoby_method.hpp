#ifndef JACOBY_METHOD_SLE_HPP
#define JACOBY_METHOD_SLE_HPP

#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"

template<typename T>

Vector<T> Jacoby_method(Matrix<T>& A, Vector<T> b, Vector<T> x, const T& tol){

    std::size_t n = A.get_rows();
    Matrix<T> C(A);

    for (std::size_t i = 0; i < n; ++i){
        for (std::size_t j = 0; j < n; ++ j){
            if (i != j)
                C[i][j] /= - C[i][i];
        }
        b[i] /= C[i][i];
        C[i][i] = 0.;
    }

    T dC = C.norm_inf();
    dC = (1 - dC) / dC;

    std::size_t iter = 0;
    Vector<T> xk(x);

    while (true) {
        ++iter;

        x = C.dot(xk) + b;  

        if ((x - xk).norm_inf() < dC * tol){
            break;
        }

        xk = x;
    }

    return x;
}

#endif