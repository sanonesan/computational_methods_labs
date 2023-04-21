#include <iostream>
#include "./include/Solver_1d_heat_eq.hpp"
#include "./include/test_equations/Test1.hpp"
#include "./include/test_equations/Test2.hpp"

int main(int args, char **argv){

    typedef long double T;
    
    Solver_1d_heat_eq<T> solver;

    Test1<T> test;
    Test2<T> test1;

    solver.solve_eq(test1);


    return 0;
}
