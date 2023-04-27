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
        T _u0 = 1.61;

        //vector = {u(x0, t), u(xL, t)}
        std::vector<std::function<T (const T x, const T t)>> _boundary_conditions;
        // 0 -- first type BC
        // 1 -- second type BC
        std::size_t _left_boundary_condition_type = 0;
        std::size_t _right_boundary_condition_type = 0;

        
        //u(x, start_time)
        std::function<T (const T x, const T t)> _initial_conditions;
        


};