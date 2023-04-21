#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_1d_heat_equation.hpp"

/*
// ...........Тест..2..из..методички.............. //
*/
template<class T>
class Test2: virtual public Class_1d_heat_equation<T>{

    public:

        Test2(){

            //material parameters
            this->_c = 5.;
            this->_rho = 1.;
            
            //Space
            this->_h = 0.1;
            this->_x0 = 0.;
            this->_xL = 1.;

            //Time
            this->_start_time = 0.;
            this->_tau = 0.005;
            this->_end_time = 0.1;
            

            //initial temperature
            this->_u0 = 30.;


            // Initial K(u, x)
            auto K = [](const T u, const T x) -> T{
                return 1. + 0.5 * x;
            };
            this->_K = K;

            // Boundary {u(x0, t), u(xL, t)}
            auto u_0_t = [this](const T x, const T t) -> T{
                return this->_u0;
            };

            auto u_L_t = [this](const T x, const T t) -> T{
                return 0.;
            };

            // 0 -- first type BC
            // 1 -- second type BC
            this->_left_boundary_condition_type = 0;
            this->_right_boundary_condition_type = 1;

            this->_boundary_conditions.push_back(u_0_t);
            this->_boundary_conditions.push_back(u_L_t);

            // Initial u(x, 0)
            auto u_x_0 = [this](const T x, const T t) -> T{
                return this->_u0 + x * (this->_xL - x);
            };

            this->_initial_conditions = u_x_0;
            
        }

};