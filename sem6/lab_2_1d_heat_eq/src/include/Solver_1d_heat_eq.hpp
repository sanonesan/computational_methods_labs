#pragma once

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>

#include <vector>
#include "Class_1d_heat_equation.hpp"
#include "./1d_heat_equation/explicit_2_layer_difference_scheme.hpp"
#include "./1d_heat_equation/implicit_2_layer_difference_scheme.hpp"
#include "./1d_heat_equation/templated_2_layer_difference_scheme.hpp"



template<class T>
class Solver_1d_heat_eq{

    private:

        void check_folder(const std::string& str){

            std::ifstream dir_stream(str.c_str());

            if (!dir_stream) {
                std::cout << "Folder created!\n" << "path:\t" << str << "\n\n";
                const int dir_err = mkdir(str.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                if (dir_err == -1)
                {
                    std::cout << ("Error creating directory!\n");
                    exit(1);
                }               
                
            }

        }

    public:

        T tol = 1e-6;
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

            this->check_folder(this->output_folder);
            std::string out_path;
            // out_path = this->output_folder + this->file_name + "explicit_1d_heat_eq_output";
            // explicit_2_layer_difference_scheme(heat_equation, out_path);
            // out_path = this->output_folder + this->file_name + "implicit_1d_heat_eq_output";
            // implicit_2_layer_difference_scheme(heat_equation, out_path);
            out_path = this->output_folder + this->file_name + "_1d_heat_eq_output";
            templated_2_layer_difference_scheme(heat_equation, (T)0.5, out_path);
            
        }


};