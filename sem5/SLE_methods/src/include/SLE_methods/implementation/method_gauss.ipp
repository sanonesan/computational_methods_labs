#include "../method_gauss.hpp"

template <typename T>
int method_gauss(Matrix<T> A, Vector<T> b, Vector<T> &solution) {
    int check = 0;

    check = straight_course(A, b);

    if (check == 0) {
        return reverse_course(A, b, solution);
    } else
        throw std::invalid_argument("\nMatrix is singular (non-invertible)\n");
};

// subfunctions1
template <typename T>
int straight_course(Matrix<T> &A, Vector<T> &b) {
    major_element(A);
    for (size_t k = 0; k < A.get_rows(); ++k) {
        if (fabs(A[k][k]) < eps) {
            std::cout << "\nMatrix is singular (non-invertible)\n";
            return 1;
        }
        remove_MatrixColumnElements_UnderLine(k, A, b);
    }
    return 0;
};

template <typename T>
int reverse_course(Matrix<T> &A, Vector<T> &b, Vector<T> &solution) {
    solution = b;
    std::size_t n = solution.size();

    solution[n - 1] /= A[n - 1][n - 1];

    for (std::size_t i = 0; i <= n - 2; ++i) {
        for (std::size_t j = n - 2 - i + 1; j < n; ++j) {
            solution[n - 2 - i] -= solution[j] * A[n - 2 - i][j];
        }
        solution[n - 2 - i] /= A[n - 2 - i][n - 2 - i];
    }

    solution.check_vector_zero();

    return 0;
};

// subfunctions2
template <typename T>
int major_element(Matrix<T> &A) {
    std::size_t i_max = 0;
    std::size_t n = A.get_rows();

    for (std::size_t k = 0; k < n; ++k) {
        i_max = k;
        for (std::size_t i = k; i < n; ++i) {
            if (fabs(fabs(A[i_max][k]) - fabs(A[i][k])) > eps) {
                if (fabs(A[i_max][k]) < fabs(A[i][k])) {
                    i_max = i;
                }
            }
        }
        swap(A[i_max], A[k]);
    }
    return 0;
};

template <typename T>
int remove_MatrixColumnElements_UnderLine(std::size_t &k, Matrix<T> &A, Vector<T> &b) {
    T tmp = 0.0;
    std::size_t n = b.size();

    for (size_t i = k + 1; i < n; ++i) {
        if (fabs(A[i][k]) > eps) {
            tmp = A[i][k] / A[k][k];

            for (size_t j = k; j < n; ++j) {
                A[i][j] -= A[k][j] * tmp;

                if (fabs(A[i][j]) < eps) {
                    A[i][j] *= (1 / A[i][j] > 0) ? (1) : (-1);
                }
            }
            b[i] -= b[k] * tmp;
            if (fabs(b[i]) < eps) {
                b[i] *= (1 / b[i] > 0) ? (1) : (-1);
            }
        } else {
            A[i][k] = 0.0;
            if (fabs(A[i][k]) < eps) {
                A[i][k] *= (1 / A[i][k] > 0) ? (1) : (-1);
            };
        }
    }

    return 0;
};
