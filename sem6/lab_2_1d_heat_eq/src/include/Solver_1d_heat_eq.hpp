#pragma once

#include<vector>


template<class T>
class Solver_1d_heat_eq{

    public:

        T tol = 1e-3;
        std::string file_name = "";
        std::string output_folder = "../output/";

        //Default constructor
        Solver_1d_heat_eq(){
            this->tol = 1e-9;
            this->output_folder = "../output/";
            this->file_name = "";
        }

        //Alt constructor
        Solver_1d_heat_eq(T tol, std::string output_folder, std::string file_name){
            this->tol = tol;
            this->output_folder = "../output/";
            this->file_name = file_name;
        }

        template<typename F>
        void solve_eq(){
            
            std::string out_path;
            out_path = this->output_folder + this->file_name + "_ode_Euler_output.csv";     
        }



};