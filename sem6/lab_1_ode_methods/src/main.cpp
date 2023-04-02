#include "./include/Solver_ODE.hpp"
#include "./include/tests/Test0.hpp"
#include "./include/tests/Test01.hpp"
#include "./include/tests/Test02.hpp"
#include "./include/tests/Test1.hpp"
#include "./include/tests/Var_4.hpp"


int main(int args, char **argv) {
    typedef long double T;
    Solver_ODE<T> solver;

    T t, t_final, tau;

    t = 0;
    t_final = 50;
    solver.tol = 1e-6;


    // APPROXIMATION ORDER TEST

    Test0<T> test0;

    solver.output_folder = "../output/test0_for_approximation_order/";

    t = 0;
    t_final = 5;
    solver.tol = 1e-5;

    //q = 0.25;
    tau = 0.01;
    solver.file_name = "test0_tau_1e-2";

    solver.solve_ode_explicit_Euler(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_2_step_symmetrical_scheme(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_RK2_fix_step(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_RK4_fix_step(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_AB4(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_Predictor_Corrector(t, t_final, tau, test0._x0, test0._ode_system);

    tau = 0.01 * 0.25;
    std::cout.precision(24);
    std::cout << tau << "\n";
    solver.file_name = "test0_tau_1e-3";

    solver.solve_ode_explicit_Euler(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_2_step_symmetrical_scheme(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_RK2_fix_step(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_RK4_fix_step(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_AB4(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_Predictor_Corrector(t, t_final, tau, test0._x0, test0._ode_system);

    tau = 0.01 * 0.25 * 0.25;
    solver.file_name = "test0_tau_1e-4";

    solver.solve_ode_explicit_Euler(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_2_step_symmetrical_scheme(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_RK2_fix_step(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_RK4_fix_step(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_AB4(t, t_final, tau, test0._x0, test0._ode_system);
    solver.solve_ode_Predictor_Corrector(t, t_final, tau, test0._x0, test0._ode_system);


    // Test01<T> test01;
    
    // tau = 0.01;
    // solver.file_name = "test0_tau_1e-2";
    // solver.solve_ode_implicit_Euler(t, t_final, tau, test01._x0, test01._ode_system);

    // tau = 0.0025;
    // solver.file_name = "test0_tau_1e-3";
    // solver.solve_ode_implicit_Euler(t, t_final, tau, test01._x0, test01._ode_system);

    // tau = 0.000625;
    // solver.file_name = "test0_tau_1e-4";
    // solver.solve_ode_implicit_Euler(t, t_final, tau, test01._x0, test01._ode_system);

    
    // APPROXIMATION ORDER TEST


    // RUNGE RULE TEST


    t = 0;
    t_final = 10;
    solver.tol = 1e-8;
    
    tau = 0.01;

    Test02<T> test02;
    solver.output_folder = "../output/test02_for_runge_rule/";
    solver.file_name = "test02_tau_1e-2";
    solver.solve_ode_RK2_vary_step(t, t_final, tau, test02._x0, test02._ode_system);
    solver.solve_ode_RK4_vary_step(t, t_final, tau, test02._x0, test02._ode_system);

    // RUNGE RULE TEST


    
    // ---------------solution_var_4--------------- //

    Var_4<T> var4;

    t = 0;
    t_final = 50;
    solver.tol = 1e-6;


    tau = 0.01;
    solver.file_name = "var4";
    solver.solve_eq_with_all_methods(t, t_final, tau, var4._x0, var4._ode_system);

    // VAR 4 PHASE PORTRAITS

    // t = 0;
    // t_final = 0.5;
    // tau = 0.01;
    // solver.tol = 1e-6;
    // solver.file_name = "var4";

    // int grid_steps = 20;


    // solver.output_folder = "../output/var4_phase_portrait/grid_step_0.1/";

    // T hi = 0.1;
    // T hj = hi;
    
    // T p_grid = 1;

    // for(std::size_t i = 0; i <= grid_steps; ++i){
    //     for(std::size_t j = 0; j <= grid_steps; ++j){
    //         var4._x0[0] = -p_grid + hi * i;
    //         var4._x0[1] = -p_grid + hj * j;
    //         solver.file_name = "var4_" + std::to_string(i) + "_" + std::to_string(j);
    //         solver.solve_ode_RK4_fix_step(t, t_final, tau, var4._x0, var4._ode_system);
    //     }
    // }

    // solver.output_folder = "../output/var4_phase_portrait/grid_step_0.05/";

    // hi = 0.05;
    // hj = hi;
    
    // p_grid = 0.5;

    // for(std::size_t i = 0; i <= grid_steps; ++i){
    //     for(std::size_t j = 0; j <= grid_steps; ++j){
    //         var4._x0[0] = -p_grid + hi * i;
    //         var4._x0[1] = -p_grid + hj * j;
    //         solver.file_name = "var4_" + std::to_string(i) + "_" + std::to_string(j);
    //         solver.solve_ode_RK4_fix_step(t, t_final, tau, var4._x0, var4._ode_system);
    //     }
    // }


    // solver.output_folder = "../output/var4_phase_portrait/grid_step_0.025/";

    // t_final = 1;
    // hi = 0.025;
    // hj = hi;
    
    // p_grid = 0.25;

    // for(std::size_t i = 0; i <= grid_steps; ++i){
    //     for(std::size_t j = 0; j <= grid_steps; ++j){
    //         var4._x0[0] = -p_grid + hi * i;
    //         var4._x0[1] = -p_grid + hj * j;
    //         solver.file_name = "var4_" + std::to_string(i) + "_" + std::to_string(j);
    //         solver.solve_ode_RK4_fix_step(t, t_final, tau, var4._x0, var4._ode_system);
    //     }
    // }

    // solver.output_folder = "../output/var4_phase_portrait/grid_step_0.0125/";

    // t_final = 1;
    // hi = 0.0125;
    // hj = hi;
    
    // p_grid = 0.125;

    // for(std::size_t i = 0; i <= grid_steps; ++i){
    //     for(std::size_t j = 0; j <= grid_steps; ++j){
    //         var4._x0[0] = -p_grid + hi * i;
    //         var4._x0[1] = -p_grid + hj * j;
    //         solver.file_name = "var4_" + std::to_string(i) + "_" + std::to_string(j);
    //         solver.solve_ode_RK4_fix_step(t, t_final, tau, var4._x0, var4._ode_system);
    //     }
    // }


    // solver.output_folder = "../output/var4_phase_portrait/grid_step_0.00625_a/";

    // t_final = 1;
    // hi = 0.01;
    // hj = hi;
    
    // for(std::size_t i = 0; i <= grid_steps; ++i){
    //     for(std::size_t j = 0; j <= grid_steps; ++j){
    //         var4._x0[0] = -0.2 + hi * i;
    //         var4._x0[1] = -0.1 + hj * j;
    //         solver.file_name = "var4_" + std::to_string(i) + "_" + std::to_string(j);
    //         solver.solve_ode_RK4_fix_step(t, t_final, tau, var4._x0, var4._ode_system);
    //     }
    // }

    // solver.output_folder = "../output/var4_phase_portrait/grid_step_0.00625_b/";

    // t_final = 1;
    // hi = 0.01;
    // hj = hi;


    // for(std::size_t i = 0; i <= grid_steps; ++i){
    //     for(std::size_t j = 0; j <= grid_steps; ++j){
    //         var4._x0[0] = -0.1 + hi * i;
    //         var4._x0[1] = -0.1 + hj * j;
    //         solver.file_name = "var4_" + std::to_string(i) + "_" + std::to_string(j);
    //         solver.solve_ode_RK4_fix_step(t, t_final, tau, var4._x0, var4._ode_system);
    //     }
    // }

    // solver.output_folder = "../output/var4_phase_portrait/grid_step_0.00625_c/";

    // t_final = 1;
    // hi = 0.01;
    // hj = hi;


    // for(std::size_t i = 0; i <= grid_steps; ++i){
    //     for(std::size_t j = 0; j <= grid_steps; ++j){
    //         var4._x0[0] = 0 + hi * i;
    //         var4._x0[1] = -0.1 + hj * j;
    //         solver.file_name = "var4_" + std::to_string(i) + "_" + std::to_string(j);
    //         solver.solve_ode_RK4_fix_step(t, t_final, tau, var4._x0, var4._ode_system);
    //     }
    // }

    // VAR 4 PHASE PORTRAITS


    // ---------------solution_var_4--------------- //

    // -------------------Test1-------------------- //

    // Test1<T> test1;

    // t = 0;
    // t_final = 10;

    // tau = 1e-2;
    // solver.tol = 1e-6;
    // solver.file_name = "test1";
    // solver.solve_eq_with_all_methods(t, t_final, tau, test1._x0, test1._ode_system);

    // -------------------Test1-------------------- //

    return 0;
}
