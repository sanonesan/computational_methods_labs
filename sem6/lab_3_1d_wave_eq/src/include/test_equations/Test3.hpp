#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_1d_wave_equation.hpp"

/**
 * Тест 3
*/
template<class T>
class Test3: virtual public Class_1d_wave_equation<T>{

    public:

        Test3(){
            
            //test name
            this->_name = std::string (__func__);

            //material parameters
            this->_a = 1.;
            
            //Space
            this->_x0 = -2.;
            this->_xL = 2.;
            this->_h = 0.01;

            //Time
            this->_start_time = 0.;
            this->_end_time = 5.;
            this->_tau = 0.005;

            // Boundary {u(x0, t), u(xL, t)}
            auto u_0_t = [this](const T x, const T t) -> T{
                return 0.;
            };

            auto u_L_t = [this](const T x, const T t) -> T{
                return 0.;
            };
            
            this->_boundary_conditions.clear();
            this->_boundary_conditions.shrink_to_fit();
            this->_boundary_conditions.push_back(u_0_t);
            this->_boundary_conditions.push_back(u_L_t);

            // Initial u(x, 0)
            auto u_x_0 = [this](const T x, const T t) -> T{
                T c1 = 1.;
                c1 /= 3;
                if(x >= -c1 && x <= c1){
                    return 1.;
                }
                else{
                    return 0.;
                }
            };

            // du/dt(x, 0)
            auto du_dt_x_0 = [this](const T x, const T t) -> T{
                return 0.;
            };

            this->_initial_conditions.clear();
            this->_initial_conditions.shrink_to_fit();
            this->_initial_conditions.push_back(u_x_0);
            this->_initial_conditions.push_back(du_dt_x_0);

        }

};