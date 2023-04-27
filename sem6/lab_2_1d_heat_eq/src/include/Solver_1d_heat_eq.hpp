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
        bool notifications = false;
        std::size_t inner_iteration_threshold = 0;
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

        template <typename sigma_type>
        void solve_eq(Class_1d_heat_equation<T>& heat_equation, const sigma_type sigma){

            this->check_folder(this->output_folder);
            std::string out_path;

            if(sigma < 0. || sigma > 1.)
                throw std::invalid_argument("sigma should be: 0 <= sigma <= 1");

            out_path = this->output_folder + this->file_name;

            if (this->notifications){
                std::cout << file_name << ": \t";
            }

            if (fabs(sigma) < 1e-16) {
                if (this->notifications){
                    std::cout << std::setw(20) << "Explicit scheme:";
                    if (heat_equation._K_type == 1){
                        std::cout << "  Warning: K = K(u, x) is nonlinear. Better chose another scheme.\t";
                    }
                }
                out_path += "_1d_heat_eq_explicit_scheme_output";             
            } else if (fabs(sigma - 1.) < 1e-16) {
                if (this->notifications){
                    std::cout << std::setw(20) << "Implicit scheme:";
                }
                out_path += "_1d_heat_eq_implicit_scheme_output";                
            } else {
                if (this->notifications){
                    std::cout << std::setw(20) << "Mixed scheme:";
                }
                out_path += "_1d_heat_eq_mixed_scheme_output";                
            }
            
            templated_2_layer_difference_scheme(heat_equation, sigma, this->inner_iteration_threshold, this->tol, out_path);
            if (this->notifications){
                std::cout << "  Done!\n";
            }
        }


};