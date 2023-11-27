#pragma once

#include<vector>
#include<cmath>
#include<functional> 

/*
// ..........Class_1d_wave_equation............... //
*/
template<class T>
class Class_2d_Poisson_equation_in_rectangle{

    public:

        std::string _name;
        
        //Space Rectangle G = [_x1_0, _x1_L1] x [_x2_0, _x2_L2]
        T _x1_0 = 0.;
        T _x1_L1 = 1.;

        T _x2_0 = 0.;
        T _x2_L2 = 1.;

        T _h1 = 0.1;
        T _h2 = 0.1;

        T _h1_2 = 0.01;
        T _h2_2 = 0.01;

        T _start_time = 0.;
        T _end_time = 5.;
        T _tau = 0.1;
        T _tau_2 = 0.01;

        // Dirichlete = 1 / Neuman = 2 boundary conditions pair = {condition type, condition function}
        std::pair<std::size_t, std::function<T (const T x1, const T x2,  const T t)>> _left_boundary_condition;
        std::pair<std::size_t, std::function<T (const T x1, const T x2,  const T t)>> _right_boundary_condition;
        std::pair<std::size_t, std::function<T (const T x1, const T x2,  const T t)>> _lower_boundary_condition;
        std::pair<std::size_t, std::function<T (const T x1, const T x2,  const T t)>> _upper_boundary_condition;
        
        
        std::function<T (const T x1, const T x2,  const T t)> _f;


        Class_2d_Poisson_equation_in_rectangle<T> DEFAULT_TEST(){

            this->_x1_0 = 0.;
            this->_x1_L1 = 1.;

            this->_x2_0 = 0.;
            this->_x2_L2 = 1.;

            this->_h1 = 0.1;
            this->_h2 = 0.1;

            this->_h1_2 = 0.01; //h_1 in power 2
            this->_h2_2 = 0.01; //h_2 in power 2

            this->_start_time = 0.;
            this->_end_time = 10.;
            this->_tau = 0.1;
            this->_tau_2 = 0.01;


            this->_left_boundary_condition.first = 1;
            this->_right_boundary_condition.first = 1;
            this->_upper_boundary_condition.first = 1;
            this->_lower_boundary_condition.first = 1;

            auto l_b = [](T x1, T x2, T t){
                return (T)1.;
            };

            auto l_b_0 = [](T x1, T x2, T t){
                return (T)1.;
            };

            this->_left_boundary_condition.second = l_b_0;
            this->_right_boundary_condition.second = l_b_0;
            this->_upper_boundary_condition.second = l_b;
            this->_lower_boundary_condition.second = l_b_0;
            
            
            auto f = [](T x1, T x2, T t){
                return (T)0.;
            };

            this->_f = f;

            return *this;           
        };

};