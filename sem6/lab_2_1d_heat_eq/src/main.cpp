#include <iostream>
#include "./include/Solver_1d_heat_eq.hpp"
#include "./include/test_equations/Test1.hpp"

int main(int args, char **argv){

    typedef double T;
    
    Solver_1d_heat_eq<T> solver;

    T t, t_final, tau, x, end_x, h;
    
    t = 0.;
    t_final = 1.;
    tau = 0.01;

    x = 0;
    end_x = 1.;
    h = 0.02;


    Test1<T> test;
    
    solver.solve_eq(t, t_final, tau, test.x0, test.xl, h, test._system);


    return 0;
}
