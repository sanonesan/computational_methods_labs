#include <iostream>

#include "./include/Solver_Fred.hpp"
#include "./include/test_equations/Test1.hpp"
// #include "./include/test_equations/Test2.hpp"

int main(int args, char **argv){

    typedef double T;

    Solver_Fred<T> solver;
    solver.notifications = true;

    Test1<T> test1;
    solver.output_folder = "../output/Test1/";
    solver.file_name = test1._name;

    solver.solve_eq_quadrature(test1);
    std::cout << "\n";

    solver.solve_eq_simple_iter(test1);
    std::cout << "\n";


    test1.DEFAULT_TEST_singular();

    solver.solve_eq_singular(test1);
    std::cout << "\n";

    return 0;
}