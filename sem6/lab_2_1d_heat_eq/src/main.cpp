#include <iostream>
#include "./include/Solver_1d_heat_eq.hpp"
#include "./include/test_equations/Test1.hpp"

int main(int args, char **argv){

    typedef double T;
    
    Solver_1d_heat_eq<T> solver;

    Test1<T> test;

    solver.solve_eq(test);


    return 0;
}
