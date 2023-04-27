#include <iostream>

#include "./include/Solver_1d_wave_eq.hpp"
#include "./include/test_equations/Test1.hpp"

int main(int args, char **argv){

    typedef long double T;
    
    Solver_1d_wave_eq<T> solver;
    solver.notifications = true;

    Test1<T> test1;
    solver.output_folder = "../output/Test1/";
    solver.file_name = test1._name;

    solver.solve_eq(test1);
    // solver.solve_eq(test1, 1);
    // solver.solve_eq(test1, 0.5);
    std::cout << "\n";

    

    return 0;
}
