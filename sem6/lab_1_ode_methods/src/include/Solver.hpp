#pragma once

#include<vector>

#include"./ode_methods/ode_Euler.hpp"
#include"./ode_methods/ode_RK4.hpp"
#include"./ode_methods/ode_Adams.hpp"

template<typename T>
class Solver{

    public:

        T tol = 0.000001;
        std::string file_name = "";
        std::string output_folder = "../output/";

        //Default constructor
        Solver(){
            this->tol = 0.000001;
            this->output_folder = "../output/";
            this->file_name = "";
        }

        //Alt constructor
        Solver(T tol, std::string output_folder, std::string file_name){
            this->tol = tol;
            this->output_folder = "../output/";
            this->file_name = file_name;
        }

        template<typename F>
        void solve_all(T t_start, T t_final, T tau, std::vector<T>&_x0, std::vector<F> &_ode_system){
            
            std::string out_path;

            // -------------------EULER-------------------- //            

            out_path = this->output_folder + this->file_name + "_ode_Euler_output.csv";            
            ode_Euler(t_start, t_final, tau, _x0, _ode_system, out_path);
            
            // -------------------EULER-------------------- //
            

            // ----------------RUNGE-KUTTA----------------- //

            out_path = this->output_folder + this->file_name + "_ode_RK4_fix_step_output.csv";
            ode_RK4_fix_step(t_start, t_final, tau, _x0, _ode_system, out_path);

            out_path = this->output_folder + this->file_name + "_ode_RK4_vary_step_output.csv";
            ode_RK4_vary_step(t_start, t_final, tau, _x0, _ode_system, this->tol, out_path);
            
            // ----------------RUNGE-KUTTA----------------- //
            

            // --------------ADAMS-BASHFORT---------------- //

            out_path = this->output_folder + this->file_name + "_ode_AB4_output.csv";
            ode_AB4(t_start, t_final, tau, _x0, _ode_system, out_path);

            // --------------ADAMS-BASHFORT---------------- //


            // ------------PREDICTOR-CORRECTOR------------- //

            out_path = this->output_folder + this->file_name + "_ode_Predictor_Corrector_output.csv";
            ode_Predictor_Corrector(t_start, t_final, tau, _x0, _ode_system, out_path);

            // ------------PREDICTOR-CORRECTOR------------- //

        }

};