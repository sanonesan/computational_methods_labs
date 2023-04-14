#pragma once

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>

#include <vector>
#include "Class_1d_heat_equation.hpp"
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

        //template<typename F>
        void solve_eq(Class_1d_heat_equation<T>& heat_equation){


            //const std::string dir("diag/curr/rq");

            std::ifstream dir_stream(this->output_folder.c_str());

            if (!dir_stream) {
                std::cout << "Folder created!\n" << "path:\t" << this->output_folder << "\n\n";
                const int dir_err = mkdir(this->output_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                if (dir_err == -1)
                {
                    std::cout << ("Error creating directory!\n");
                    exit(1);
                }               
                
            }

            std::string out_path;
            out_path = this->output_folder + this->file_name + "_1d_heat_eq_output";
            scheme_1(heat_equation, out_path);
        }


};