#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_1d_heat_equation.hpp"

/**
 * Линейная среда
 * Постоянная температура на границах
 * K(u, x) = const
 * 
*/
template<class T>
class Test1: virtual public Class_1d_heat_equation<T>{

    public:

        Test1(){

            //material parameters
            this->_c = 1.;
            this->_rho = 1.;
            
            //Space
            this->_x0 = 0.;
            this->_xL = 1.;
            this->_h = 0.1;

            //Time
            this->_start_time = 0.;
            this->_end_time = 0.1;
            this->_tau = 0.005;

            //initial temperature
            this->_u0 = 10.;


            // Initial K(u, x)
            auto K = [](const T u, const T x) -> T{
                return 1.;
            };
            this->_K = K;
            this->_K_type = 0;


            // Boundary {u(x0, t), u(xL, t)}
            auto u_0_t = [this](const T x, const T t) -> T{
                return this->_u0;
            };

            auto u_L_t = [this](const T x, const T t) -> T{
                return this->_u0;
            };

            this->_left_boundary_condition_type = 0;
            this->_right_boundary_condition_type = 0;

            this->_boundary_conditions.push_back(u_0_t);
            this->_boundary_conditions.push_back(u_L_t);

            // Initial u(x, 0)
            auto u_x_0 = [this](const T x, const T t) -> T{
                return this->_u0 + x * (this->_xL - x);
            };

            this->_initial_conditions = u_x_0;
            
        }

};