#include <iostream>
#include "./include/Solver_1d_heat_eq.hpp"
#include "./include/test_equations/Test1.hpp"
#include "./include/test_equations/Test2.hpp"
#include "./include/test_equations/Test3.hpp"

int main(int args, char **argv){

    typedef long double T;
    
    Solver_1d_heat_eq<T> solver;

    Test1<T> test1;
    Test2<T> test2;
    Test3<T> test3;

    solver.solve_eq(test2);


    return 0;
}
