#ifndef QR_METHOD_HPP
#define QR_METHOD_HPP

#include <iostream>

#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"

template <typename T>
int QR_method(
    Matrix<T> A,
    Vector<T> b,
    Matrix<T> &Q,
    Matrix<T> &R,
    Vector<T> &solution);

template <typename T>
int Rb_sys_matrix(
    const Matrix<T> &A,
    Vector<T> &b,
    Matrix<T> &Q,
    Matrix<T> &R);

template <typename T>
int coefs(const std::size_t &k, const std::size_t &l, T &c, T &s, Matrix<T> &A);

#include "./implementation/QR_method.ipp"

#endif
