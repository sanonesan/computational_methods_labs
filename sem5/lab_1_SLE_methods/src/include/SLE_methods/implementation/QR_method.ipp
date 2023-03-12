#include <cmath>
#include <iostream>

#include "../QR_method.hpp"
#include "../method_gauss.hpp"

#define eps 1e-16

template <typename T>
int QR_method(
    Matrix<T> A,
    Vector<T> b,
    Matrix<T> &Q,
    Matrix<T> &R,
    Vector<T> &solution) {
    Rb_sys_matrix(A, b, Q, R);

    for (std::size_t k = 0; k < R.get_rows(); ++k)
        if (fabs(R[k][k]) < eps) {
            std::cout << "\nMatrix is singular (non-invertible)\n";
            return 1;
        }

    reverse_course(R, b, solution);

    return 0;
}

template <typename T>
int Rb_sys_matrix(
    Matrix<T> &A,
    Vector<T> &b_vec,
    Matrix<T> &Q,
    Matrix<T> &R) {
    std::size_t n = A.get_rows();

    T c = 0., s = 0.;

    Q.make_matrix_identity(n);
    R = A;

    T a = 0., b = 0.;

    for (std::size_t i = 0; i < n - 1; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            coefs(i, j, c, s, R);

            for (std::size_t k = 0; k < n; ++k) {
                a = R[i][k];
                b = R[j][k];
                R[i][k] = (c * a + s * b);
                R[j][k] = (-s * a + c * b);

                a = Q[i][k];
                b = Q[j][k];
                Q[i][k] = (c * a + s * b);
                Q[j][k] = (-s * a + c * b);
            }

            a = b_vec[i];
            b = b_vec[j];
            b_vec[i] = (c * a + s * b);
            b_vec[j] = (-s * a + c * b);
        }
    }

    Q.check_matrix_zero();
    R.check_matrix_zero();
    b_vec.check_vector_zero();

    return 0;
}

template <typename T>
int coefs(
    const std::size_t &k,
    const std::size_t &l,
    T &c,
    T &s,
    Matrix<T> &A) {
    std::size_t n = A.get_rows();

    if (k < n && l < n) {
        T temp = sqrt(pow(A[k][k], 2) + pow(A[l][k], 2));

        c = A[k][k] / temp;

        s = A[l][k] / temp;

        if (fabs(c) < eps)
            c = (1 / c > 0) ? c : c * (-1);
        if (fabs(s) < eps)
            s = (1 / s > 0) ? s : s * (-1);
    }

    return 0;
}
