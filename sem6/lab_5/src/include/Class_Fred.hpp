#pragma once

#include <cmath>
#include <functional>
#include <vector>

/*
// ..........Class_1d_wave_equation............... //
*/
template <class T>
class Class_Fred {
    public:
        std::string _name;

        /**
         * For usual
         */

        T _x_0 = 0.;
        T _x_L = 1.;

        T _h = 0.01;

        T _lambda = 1;

        std::size_t _Invertible = 0;
        std::size_t _Singular = 0;
        std::function<T(const T _x, const T _s)> _K;
        std::function<T(const T _x)> _f;
        /**
         * For singular
         */

        T delta_l = 1.;
        std::vector<std::pair<T, T>> k_i;
        std::vector<std::pair<T, T>> c_i;
        std::vector<std::pair<T, T>> n_i;
        std::function<std::pair<T, T>(std::pair<T, T>, std::pair<T, T>)> _Q;
        std::function<T(std::pair<T, T>)> _f_singular;

        Class_Fred<T> DEFAULT_TEST(T x_0 = 0)  {
            this->_x_0 = x_0;
            // this->_x_L = M_PI;
            this->_x_L = 1.;
            this->_h = 0.1;
            this->_lambda = 1. / 2;
            // this->_lambda = 1. / 2 / M_PI;
            this->_Invertible = 0;
            this->_Singular = 0;

            auto K = [](T _x, T _s) -> T {
                return (1 - _x * cos(_x * _s));
            };

            auto f = [](T _x) -> T {
                return (1 + sin(_x)) / 2;
            };

            // auto f = [](T x) -> T {
            //     return  x*x + sqrt(x);
            // };

            // auto f = [](T _x) -> T {
            //     return  cos(_x);
            // };

            this->_K = K;
            this->_f = f;

            return *this;
        };

        Class_Fred<T> DEFAULT_TEST_1(T x_0 = 0) {
            this->_x_0 = x_0;
            // this->_x_L = M_PI;
            this->_x_L = 1.;
            this->_h = 0.1;
            this->_lambda = 1. / 2;
            // this->_lambda = 1. / 2 / M_PI;
            this->_Invertible = 0;
            this->_Singular = 0;

            auto K = [](T _x, T _s) -> T {
                return (1 - _x * cos(_x * _s));
            };

            // auto K = [](T _x, T _s) -> T {
            //     return (_s * sin(_x));
            // };

            // auto f = [](T _x) -> T {
            //     return  (1 + sin(_x)) / 2;
            // };

            auto f = [](T x) -> T {
                return x * x + sqrt(x);
            };

            // auto f = [](T _x) -> T {
            //     return  cos(_x);
            // };

            this->_K = K;
            this->_f = f;

            return *this;
        };

        Class_Fred<T> DEFAULT_TEST_singular(int N = 10) {

            this->delta_l = 2 * M_PI / N;

            this->_Singular = 1;
            std::pair<T, T> tmp = {0., 0.};

            for (std::size_t i = 1; i <= N; ++i) {
                tmp.first = cos(2 * M_PI * (i - 0.5) / N);
                tmp.second = sin(2 * M_PI * (i - 0.5) / N);
                this->k_i.push_back(tmp);
                this->n_i.push_back(tmp);

                tmp.first = cos(2 * M_PI * (i - 1) / N);
                tmp.second = sin(2 * M_PI * (i - 1) / N);
                this->c_i.push_back(tmp);
            }

            auto Q = [](std::pair<T, T> k_i, std::pair<T, T> c_i) -> std::pair<T, T> {
                std::pair<T, T> tmp = {0., 0.};

                tmp.first = (-(k_i.second - c_i.second)) / 2 / M_PI / (pow((k_i.first - c_i.first), 2) + pow((k_i.second - c_i.second), 2));

                tmp.second = (k_i.first - c_i.first) / 2 / M_PI / (pow((k_i.first - c_i.first), 2) + pow((k_i.second - c_i.second), 2));

                return tmp;
            };

            this->_Q = Q;

            this->_f_singular = [](std::pair<T, T> r) -> T {
                T phi = 0;

                if (r.first >= 0 && r.second >= 0) {
                    phi = atan(r.second / r.first);
                } else if (r.first <= 0 && r.second >= 0) {
                    phi = M_PI + atan(r.second / r.first);
                } else if (r.first <= 0 && r.second <= 0) {
                    phi = M_PI + atan(r.second / r.first);
                } else if (r.first >= 0 && r.second <= 0) {
                    phi = 2 * M_PI + atan(r.second / r.first);
                }

                return sin(3 * phi);
            };

            return *this;
        };
    };