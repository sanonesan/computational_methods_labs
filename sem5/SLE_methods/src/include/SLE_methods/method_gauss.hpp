#ifndef METHOD_GAUSS_HPP
#define METHOD_GAUSS_HPP

#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"

template <typename T>
int method_gauss(Matrix<T> A, Vector<T> b, Vector<T> &solution);

template <typename T>
int straight_course(Matrix<T> &A, Vector<T> &b);

template <typename T>
int reverse_course(Matrix<T> &A, Vector<T> &b, Vector<T> &solution);

template <typename T>
int major_element(Matrix<T> &A);

template <typename T>
int remove_MatrixColumnElements_UnderLine(std::size_t &k, Matrix<T> &A, Vector<T> &b);

#include "./implementation/method_gauss.ipp" 

#endif
