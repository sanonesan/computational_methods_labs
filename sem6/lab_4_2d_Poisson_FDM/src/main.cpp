#include <iostream>

#include "./include/Solver_2d_Poisson_eq.hpp"
#include "./include/test_equations/Test1.hpp"
// #include "./include/test_equations/Test2.hpp"

int main(int args, char **argv){

    typedef double T;
    
    Solver_2d_Poisson_eq<T> solver;
    solver.notifications = true;

    Test1<T> test1;
    solver.output_folder = "../output/Test1/";
    solver.file_name = test1._name;

    solver.solve_eq(test1);
    std::cout << "\n";

    return 0;
}
