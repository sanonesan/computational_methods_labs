#ifndef SIMPLE_ITER_METHOD_SLE_HPP
#define SIMPLE_ITER_METHOD_SLE_HPP

#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"

template<typename T>

Vector<T> simple_iter_method(Matrix<T>& A, Vector<T> b, Vector<T> x, const T& tau, const std::size_t& M, const T& tol){

    std::size_t n = A.get_rows();
    Matrix<T> C(A);

    C *= - tau;

    for (std::size_t i = 0; i < n; ++i){
        C[i][i] += 1.;
        b[i] *= tau;
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

        if (iter == M && M != 0)
            break;

        xk = x;
    }
    std::cout << iter << "\n";
    return x;
}

#endif