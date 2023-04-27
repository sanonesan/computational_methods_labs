#pragma once

#include<vector>
#include<cmath>
#include<functional> 

/*
// ..........Class_1d_heat_equation............... //
*/
template<class T>
class Class_1d_heat_equation{

    public:

        std::string _name;

        //material parameters
        T _c = 1.;
        T _rho = 1.;
        std::function<T (const T u, const T x)> _K;
        // 0 -- linear K (K = const || K = K(x))
        // 1 -- nonlinear K (K = K(u) || K = K(u, x))
        std::size_t _K_type = 0;

        //Space
        T _x0 = 0.;
        T _xL = 1.;
        T _h = 0.1;

        //Time
        T _start_time = 0.;
        T _end_time = 0.1;
        T _tau = 0.005;

        //initial temperature
        T _u0 = 0.;

        //vector = {u(x0, t), u(xL, t)}
        std::vector<std::function<T (const T x, const T t)>> _boundary_conditions;
        // 0 -- first type BC
        // 1 -- second type BC
        std::size_t _left_boundary_condition_type = 0;
        std::size_t _right_boundary_condition_type = 0;

        
        //u(x, start_time)
        std::function<T (const T x, const T t)> _initial_conditions;
        

        /**
         * Custom constructor
         * 
         * Default test:
         * c = 1.;
         * rho = 1.;
         * 
         * x0 = 0.;
         * xL = 1.;
         * _h = 0.1;
         * 
         * start_time = 0.;
         * end_time = 0.1;
         * tau = 0.005;
         * 
         * u0 = 0.;
         * 
         * K = 1. = const;
         * K_type = 0;
         * 
         * Boundary contitions:
         * u_0_t = u0;
         * u_L_t = u0;         * 
         * left_boundary_condition_type = 0;
         * right_boundary_condition_type = 0;
         * 
         * Initial contitions:
         * u_x_0 = u0 + x * (xL - x);
         * 
        */
        Class_1d_heat_equation<T> DEFAULT_TEST(){

            //test name
            this->_name = std::string (__func__);


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
            this->_u0 = 0.;


            // Initial K(u, x)
            auto K = [](const T u, const T x) -> T{
                return 1.;
            };
            this->_K = K;
            // 0 -- linear K (K = const || K = K(x))
            // 1 -- nonlinear K (K = K(u) || K = K(u, x))
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
            
        };

};