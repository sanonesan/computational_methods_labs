#include"./include/Solver.hpp"

#include"./include/tests/Var_4.hpp"
#include"./include/tests/Test1.hpp"



int main(int args, char **argv){

	typedef double T;
	Solver<T> solver;

	T t, t_final, tau;	

	// ---------------solution_var_4--------------- //
	
	Var_4<T> var4;	
	
	t = 0;
	t_final = 50;
	tau = 0.001;

	solver.tol = 0.001;
	solver.file_name = "var4";
	solver.solve_all(t, t_final, tau, var4._x0, var4._ode_system);	

	// ---------------solution_var_4--------------- //


	// -------------------Test1-------------------- //

	Test1<T> test1;

	t = 0;
	t_final = 20;
	tau = 0.01;
	
	solver.tol = 0.01;
	solver.file_name = "test1";
	solver.solve_all(t, t_final, tau, test1._x0, test1._ode_system);

	// -------------------Test1-------------------- //

    return 0;
}

