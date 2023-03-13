#pragma once

#include <vector>
#include <cmath>
#include <string>


#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"

template<typename T>
using matrix = std::vector<std::vector<T>>;

// умножение матрицы и столбца
template<typename T>
std::vector<T> mult(const matrix<T>& A1, const std::vector<T>& A2) {
	std::size_t n = A2.size();
	std::vector<T> res;
	res.reserve(n); 
	for (std::size_t i = 0; i < n; ++i) {
		T c = 0;
		for (std::size_t j = 0; j < n; ++j) {
			c += A1[i][j] * A2[j];
		}
		res.push_back(c);

	}
	return res;
}

// Обращение матрицы 2x2 (просто ввел формулы из Wolfram Mathematica)
template<typename T>
matrix<T> inv(const matrix<T>& m) {
	T det = m[0][0] * m[1][1] - m[1][0] * m[0][1];
	return matrix<T> {
		{  m[1][1] / det, -m[0][1] / det },
		{ -m[1][0] / det,  m[0][0] / det }
	};
}

// Численная производная функции
template<typename T, typename F>
auto num_deriv(F& foo, T _eps) {
	return [&foo, _eps](T x) {
		return (foo(x + _eps) - foo(x)) / _eps;
	};
}

template<typename T, typename Func, typename J>
std::vector<T>  newton_sys_method1(Func foo1, Func foo2, J jacobi, const std::vector<T>& x0, std::size_t &iter_counter)
{
	std::vector<T> xk = x0;
	T coef;

	iter_counter = 0;

	do {

		matrix<T> jacobiInv = inv(jacobi(xk[0], xk[1]));

		std::vector<T> F{ foo1(xk[0], xk[1]), foo2(xk[0], xk[1]) };

		std::vector<T> jacobiInvMultF = mult(jacobiInv, F);

		std::vector<T> newXk;

		newXk.reserve(x0.size());
		for (std::size_t i = 0; i < x0.size(); ++i) {
			newXk.push_back(xk[i] - jacobiInvMultF[i]);
		}

		coef = norm(jacobiInvMultF);

		xk = newXk;

		++iter_counter;

	} while (coef > eps);

	//cout << "Количество итераций: " << iter_counter << "\n";

	return xk;
}


/*
    Невный метод Эйлера
*/
template <typename T, typename F>
void ode_implicit_Euler(T start_time, T end_time, T tau, std::vector<T> x, const std::vector<F> &func, const std::string &out_path, const T &tol) {
    std::ofstream fout(out_path);
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
    fout << "time";
    for (std::size_t i = 0; i < func.size(); ++i) {
        fout << ",u" << i;
    }
    fout << "\n";
    fout << start_time;
    for (std::size_t i = 0; i < func.size(); ++i) {
        fout << "," << x[i];
    }
    fout << "\n";

    std::vector<T> tmp(x);

    while (start_time <= end_time) {
        fout << start_time + tau;

        for (std::size_t i = 0; i < x.size(); ++i) {

            /**
             *
             *
             * ????????????????
             *
             *
             *
             */

            fout << "," << x[i];
        }
        fout << "\n";
        tmp.assign(x.begin(), x.end());

        start_time += tau;
    }

    return;
}
