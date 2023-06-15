#include <iostream>

#include "./include/Solver_Fred.hpp"
#include "./include/test_equations/Test1.hpp"
#include "./include/test_equations/Test2.hpp"
#include "./include/test_equations/Test3.hpp"
// #include "./include/test_equations/Test2.hpp"

int main(int args, char **argv){

    typedef double T;

    Solver_Fred<T> solver;
    solver.notifications = true;

    Test1<T> test1;
    solver.output_folder = "../output/Test1/";
    solver.file_name = test1._name;
    solver.solve_eq_quadrature(test1);
    solver.solve_eq_simple_iter(test1);

    test1.set_01();
    solver.file_name = test1._name;
    solver.solve_eq_quadrature(test1);
    solver.solve_eq_simple_iter(test1);
    std::cout << "\n";


    std::cout << "\n";

    Test2<T> test2;
    solver.output_folder = "../output/Test2/";
    solver.file_name = test2._name;
    solver.solve_eq_quadrature(test2);
    solver.solve_eq_simple_iter(test2);

    test2.set_01();
    solver.file_name = test2._name;
    solver.solve_eq_quadrature(test2);
    solver.solve_eq_simple_iter(test2);
    std::cout << "\n";

    Test3<T> test3(100);
    solver.output_folder = "../output/Test3/";
    solver.file_name = test3._name;
    solver.solve_eq_singular(test3);
    std::cout << "\n";

    return 0;
}