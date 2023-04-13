#include <iostream>
#include "./include/Solver_1d_heat_eq.hpp"
#include "./include/test_equations/Test1.hpp"

int main(int args, char **argv){

    typedef double T;
    
    Solver_1d_heat_eq<T> solver;

    T t, t_final, tau, x, end_x, h;
    
    t = 0.;
    t_final = 5.;
    tau = 0.1;

    x = 0.;
    end_x = 5.;
    h = 0.1;

    Test1<T> test;
    
    solver.solve_eq(t, t_final, tau, x, end_x, h, test._system);


    return 0;
}
