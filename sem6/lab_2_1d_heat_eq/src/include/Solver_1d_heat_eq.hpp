#pragma once

#include <vector>
#include "./1d_heat_equation/scheme_1.hpp"

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
        void solve_eq(T start_time, T end_time, T tau, T start_x, T end_x, T h, const std::vector<F> &func){

            std::string out_path;
            out_path = this->output_folder + this->file_name + "_1d_heat_eq_output.csv";
            scheme1(start_time, end_time, tau, start_x, end_x, h, func, out_path);
        }


};