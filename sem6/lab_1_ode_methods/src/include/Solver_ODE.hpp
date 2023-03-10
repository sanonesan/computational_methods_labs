#pragma once

#include<vector>

#include"./ode_methods/ode_Euler.hpp"
#include"./ode_methods/ode_RK4.hpp"
#include"./ode_methods/ode_Adams.hpp"



/*
    Solver_ODE --- класс для решения ODE методами:
        - явный метод Эйлера;
        - метод Рунге-Кутты 4 порядка с фиксированным шагом;
        - метод Рунге-Кутты 4 порядка с изменяющимя шагом;
        - метод Адамса-Башфорта (первые 3 шага --- Рунге-Кутта 4 порядка);
        - метод прогноз-коррекция (первые 3 шага --- Рунге-Кутта 4 порядка).
*/
template<class T>
class Solver_ODE{

    public:

        T tol = 1e-9;
        std::string file_name = "";
        std::string output_folder = "../output/";

        //Default constructor
        Solver_ODE(){
            this->tol = 1e-9;
            this->output_folder = "../output/";
            this->file_name = "";
        }

        //Alt constructor
        Solver_ODE(T tol, std::string output_folder, std::string file_name){
            this->tol = tol;
            this->output_folder = "../output/";
            this->file_name = file_name;
        }


        /*
            Явный метод Эйлера
        */
        template<typename F>
        void solve_ode_Euler(T t_start, T t_final, T tau, std::vector<T>&_x0, std::vector<F> &_ode_system){

            std::string out_path;            
            out_path = this->output_folder + this->file_name + "_ode_Euler_output.csv";            
            ode_Euler(t_start, t_final, tau, _x0, _ode_system, out_path);
            
        
        }


        /*
            Метод Рунге-Кутты 4 порядка с фиксированным шагом
        */
        template<typename F>
        void solve_ode_RK4_fix_step(T t_start, T t_final, T tau, std::vector<T>&_x0, std::vector<F> &_ode_system){

            std::string out_path;            
            out_path = this->output_folder + this->file_name + "_ode_RK4_fix_step_output.csv";
            ode_RK4_fix_step(t_start, t_final, tau, _x0, _ode_system, out_path);

        }


        /*
            Метод Рунге-Кутты 4 порядка с изменяющимя шагом
        */
        template<typename F>
        void solve_ode_RK4_vary_step(T t_start, T t_final, T tau, std::vector<T>&_x0, std::vector<F> &_ode_system){

            std::string out_path;            
            out_path = this->output_folder + this->file_name + "_ode_RK4_vary_step_output.csv";
            ode_RK4_vary_step(t_start, t_final, tau, _x0, _ode_system, this->tol, out_path);

        }


        /*
            Метод Адамса-Башфорта (первые 3 шага --- Рунге-Кутта 4 порядка)
        */
        template<typename F>
        void solve_ode_AB4(T t_start, T t_final, T tau, std::vector<T>&_x0, std::vector<F> &_ode_system){

            std::string out_path;            
            out_path = this->output_folder + this->file_name + "_ode_AB4_output.csv";
            ode_AB4(t_start, t_final, tau, _x0, _ode_system, out_path);

            
        }


        /*
            Метод прогноз-коррекция (первые 3 шага --- Рунге-Кутта 4 порядка)
        */
        template<typename F>
        void solve_ode_Predictor_Corrector(T t_start, T t_final, T tau, std::vector<T>&_x0, std::vector<F> &_ode_system){
            
            std::string out_path;            
            out_path = this->output_folder + this->file_name + "_ode_Predictor_Corrector_output.csv";
            ode_Predictor_Corrector(t_start, t_final, tau, _x0, _ode_system, out_path);

        }


        /*
            Решение ODE всеми методами:
                - явный метод Эйлера;
                - метод Рунге-Кутты 4 порядка с фиксированным шагом;
                - метод Рунге-Кутты 4 порядка с изменяющимя шагом;
                - метод Адамса-Башфорта (первые 3 шага --- Рунге-Кутта 4 порядка);
                - метод прогноз-коррекция (первые 3 шага --- Рунге-Кутта 4 порядка).
        */
        template<typename F>
        void solve_eq_with_all_methods(T t_start, T t_final, T tau, std::vector<T>&_x0, std::vector<F> &_ode_system){
            
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
