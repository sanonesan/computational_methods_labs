#ifndef NUMERIACAL_DERIVATIVE_FOR_ODE_SYS_HPP
#define NUMERIACAL_DERIVATIVE_FOR_ODE_SYS_HPP

#include <functional>
#include <vector>


#define eps 1e-16


/*
        Численно найденная матрица Якоби в точке
*/
template <typename T>
T numerical_derivative_for_ode_sys(std::vector<std::function<T(const std::vector<T> &, const T)>> func, std::vector<T> &x, T time, std::size_t i) {
    std::vector<T> tmp(x);
    T _eps = 1e-9;
    tmp[i] = x[i] + _eps;
    return (func[i](tmp, time) - func[i](x, time)) / _eps;
}


// /*
//         Численно найденная матрица Якоби в точке
//         + запоминание значений в точке х
// */
// template <typename T>
// Matrix<T> numerical_Jacobian_for_ode_sys_memorize(std::vector<std::function<T(const std::vector<T> &, const T)>> func, std::vector<T> &x, T time, std::vector<T> &foo) {
//     Matrix<T> Jacobian(func.size(), x.size());
//     std::vector<T> tmp(x);

//     T _eps = 1e-9;

//     for (std::size_t i = 0; i < Jacobian.get_rows(); ++i) {
//         for (std::size_t j = 0; j < Jacobian.get_cols(); ++j) {
//             tmp[j] = x[j] + _eps;
//             foo[i] = func[i](x, time);
//             Jacobian[i][j] = (func[i](tmp, time) - foo[i]) / _eps;
//             tmp[j] = x[j];
//         }
//     }

//     return Jacobian;
// }

#endif