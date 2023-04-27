#include <iostream>

#include "./include/Solver_1d_heat_eq.hpp"
#include "./include/test_equations/Test1.hpp"
#include "./include/test_equations/Test2.hpp"
#include "./include/test_equations/Test3.hpp"
#include "./include/test_equations/Test4.hpp"
#include "./include/test_equations/Test5.hpp"

int main(int args, char **argv){

    typedef double T;
    
    Solver_1d_heat_eq<T> solver;

    // Test1<T> test1;
    // solver.output_folder = "../output/Test1/";
    // solver.file_name = "Test1";

    // solver.solve_eq(test1, 0);
    // solver.solve_eq(test1, 1);
    // solver.solve_eq(test1, 0.5);

    // solver.output_folder = "../output/Test1/";
    // solver.file_name = "Test1_not_mono";
    
    // test1._tau = 0.01;
    // test1._h;
    // solver.solve_eq(test1, 0);
    // solver.solve_eq(test1, 1);
    // solver.solve_eq(test1, 0.5);


    // Test2<T> test2;
    // solver.output_folder = "../output/Test2/";
    // solver.file_name = "Test2";

    // solver.solve_eq(test2, 0);
    // solver.solve_eq(test2, 1);
    // solver.solve_eq(test2, 0.5);



    Test3<T> test3;
    solver.output_folder = "../output/Test3/";
    solver.file_name = "Test3";
    solver.notifications = true;
    test3._K_type = 1;
    solver.solve_eq(test3, 0);
    solver.solve_eq(test3, 1);
    solver.solve_eq(test3, 0.5);



    // Test4<T> test4;
    // solver.output_folder = "../output/Test4/";
    // solver.file_name = "Test4";

    // solver.solve_eq(test4, 0);
    // solver.solve_eq(test4, 1);
    // solver.solve_eq(test4, 0.5);


    // Test5<T> test5;
    // solver.output_folder = "../output/Test5/";
    // solver.file_name = "Test5";

    // solver.solve_eq(test5, 0);
    // solver.solve_eq(test5, 1);
    // solver.solve_eq(test5, 0.5);

    return 0;
}
