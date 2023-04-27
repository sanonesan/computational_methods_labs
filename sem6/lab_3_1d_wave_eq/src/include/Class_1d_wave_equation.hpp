#pragma once

#include<vector>
#include<cmath>
#include<functional> 

/*
// ..........Class_1d_wave_equation............... //
*/
template<class T>
class Class_1d_wave_equation{

    public:

        std::string _name;

        //Velocity of small disturbances.
        T _a = 1.;
        
        //Space
        T _x0 = 0.;
        T _xL = 1.;
        T _h = 0.1;

        //Time
        T _start_time = 0.;
        T _end_time = 0.1;
        T _tau = 0.005;

        //vector = {u(x0, t), u(xL, t)}
        std::vector<std::function<T (const T x, const T t)>> _boundary_conditions;
                
        //vector = {u(x, start_time), du/dt(x, start_time)}
        std::vector<std::function<T (const T x, const T t)>> _initial_conditions;
        

        /**
         * Custom constructor
         * 
         * Default test:
         * a = 1.;
         * 
         * x0 = 0.;
         * xL = 1.;
         * _h = 0.1;
         * 
         * start_time = 0.;
         * end_time = 0.1;
         * tau = 0.005;
         * 
         * Boundary contitions:
         * u_0_t = 0.;
         * u_L_t = 0.;         
         * 
         * Initial contitions:
         * u_x_0 = sin( M_PI * x);
         * du_dt_x_0 = 0.;
         * 
        */
        Class_1d_wave_equation<T> DEFAULT_TEST(){

            //material parameters
            this->_a = 1.;
            
            //Space
            this->_x0 = 0.;
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
                return sin(M_PI * x);
            };

            auto ddu_dxx_x_0 = [this](const T x, const T t) -> T{
                return -sin(M_PI * x);
            };

            // du/dt(x, 0)
            auto du_dt_x_0 = [this](const T x, const T t) -> T{
                return 0.;
            };

            this->_initial_conditions.clear();
            this->_initial_conditions.shrink_to_fit();
            this->_initial_conditions.push_back(u_x_0);
            this->_initial_conditions.push_back(du_dt_x_0);
            this->_initial_conditions.push_back(ddu_dxx_x_0);

            return *this;           
        };

};