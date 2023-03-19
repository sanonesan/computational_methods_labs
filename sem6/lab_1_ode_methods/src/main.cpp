#include "./include/Solver_ODE.hpp"
#include "./include/tests/Test0.hpp"
#include "./include/tests/Test1.hpp"
#include "./include/tests/Var_4.hpp"




int main(int args, char **argv) {
    typedef double T;
    Solver_ODE<T> solver;

    T t, t_final, tau;

    
    // ---------------solution_var_4--------------- //

    Var_4<T> var4;

    t = 0;
    t_final = 50;
    solver.tol = 1e-6;


    tau = 0.01;
    solver.file_name = "var4";
    solver.solve_eq_with_all_methods(t, t_final, tau, var4._x0, var4._ode_system);

    // ---------------solution_var_4--------------- //

    // -------------------Test1-------------------- //

    Test1<T> test1;

    t = 0;
    t_final = 10;

    tau = 0.01;
    solver.tol = 1e-6;
    solver.file_name = "test1";
    solver.solve_eq_with_all_methods(t, t_final, tau, test1._x0, test1._ode_system);

    // -------------------Test1-------------------- //

    return 0;
}
