#include <iostream>
#include <chrono>

#include "./include/Solver_1d_heat_eq.hpp"
#include "./include/test_equations/Test1.hpp"
#include "./include/test_equations/Test2.hpp"
#include "./include/test_equations/Test3.hpp"
#include "./include/test_equations/Test4.hpp"
#include "./include/test_equations/Test5.hpp"

int main(int args, char **argv){

    // auto begin = std::chrono::steady_clock::now();

    typedef double T;
    
    Solver_1d_heat_eq<T> solver;
    solver.notifications = true;

    Test1<T> test1;
    solver.output_folder = "../output/Test1/";
    solver.file_name = test1._name;
    
    solver.solve_eq(test1, 0);
    solver.solve_eq(test1, 1);
    solver.solve_eq(test1, 0.5);
    std::cout << "\n";

    solver.output_folder = "../output/Test1/";
    solver.file_name = test1._name + "_not_mono";
    
    test1._tau = 0.01;
    test1._h;
    test1._end_time = 0.1;
    test1._K = [](T u, T x) -> T{
        return 100.;
    };
    solver.solve_eq(test1, 0);
    solver.solve_eq(test1, 1);
    solver.solve_eq(test1, 0.5);
    std::cout << "\n";


    // Test2<T> test2;
    // solver.output_folder = "../output/Test2/";
    // solver.file_name = test2._name;

    // solver.solve_eq(test2, 0);
    // solver.solve_eq(test2, 1);
    // solver.solve_eq(test2, 0.5);
    // std::cout << "\n";


    // Test3<T> test3;
    // solver.output_folder = "../output/Test3/";
    // solver.file_name = test3._name;
    // solver.inner_iteration_threshold = 100;
    // solver.solve_eq(test3, 0);
    // solver.solve_eq(test3, 1);
    // solver.solve_eq(test3, 0.5);
    // std::cout << "\n";
    // solver.inner_iteration_threshold = 0;
    

    // Test4<T> test4;
    // solver.output_folder = "../output/Test4/";
    // solver.file_name = test4._name;

    // solver.solve_eq(test4, 0);
    // solver.solve_eq(test4, 1);
    // solver.solve_eq(test4, 0.5);
    // std::cout << "\n";


    // Test5<T> test5;
    // solver.output_folder = "../output/Test5/";
    // solver.file_name = test5._name;

    // solver.solve_eq(test5, 0);
    // solver.solve_eq(test5, 1);
    // solver.solve_eq(test5, 0.5);
    // std::cout << "\n";

    // auto end = std::chrono::steady_clock::now();
    // auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    
    // std::cout << "\nTIME: " << elapsed_ms.count() << std::endl;
    return 0;
}
