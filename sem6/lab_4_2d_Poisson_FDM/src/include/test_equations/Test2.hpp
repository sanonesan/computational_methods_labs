#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_2d_Poisson_equation_in_rectangle.hpp"

/**
 * Тест 1
*/
template<class T>
class Test2: virtual public Class_2d_Poisson_equation_in_rectangle<T>{

    public:

        Test2(){

            this->_name = std::string (__func__);


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

            this->_upper_boundary_condition.first = 2;
            this->_lower_boundary_condition.first = 2;

            auto l_b = [](T x1, T x2, T t){
                return 1 + x2;

            };

            auto r_b = [](T x1, T x2, T t){
                return 1 + x2;
            };

            auto low_b = [](T x1, T x2, T t) -> T {
                return (T)1.;
            };

            auto up_b = [](T x1, T x2, T t) -> T {
                return (T)-1.;
            };

            this->_left_boundary_condition.second = l_b;
            this->_right_boundary_condition.second = r_b;
            this->_upper_boundary_condition.second = up_b;
            this->_lower_boundary_condition.second = low_b;
            
            
            auto f = [](T x1, T x2, T t){
                return (T)0.;
            };

            this->_f = f;

            // return *this;           
        };

        

};