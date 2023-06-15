#pragma once

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>

#include <vector>
#include "Class_Fred.hpp"
#include "./Fred/Fred_init.hpp"

template<class T>
class Solver_Fred{

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

        T tol = 1e-9;
        bool notifications = false;
        std::string file_name = "";
        std::string output_folder = "../output/";

        //Default constructor
        Solver_Fred(){
            this->tol = 1e-9;
            this->output_folder = "../output/";
            this->file_name = "";
        }

        //Alt constructor
        Solver_Fred(T tol, std::string output_folder, std::string file_name){
            this->tol = tol;
            this->output_folder = "../output/";
            this->file_name = file_name;
        }

        void solve_eq_quadrature(Class_Fred<T>& Fred_equation){

            this->check_folder(this->output_folder);
            std::string out_path;

            out_path = this->output_folder + this->file_name;

            if (this->notifications){
                std::cout << file_name << ": \t";
            }

            out_path += "_Fred_quadrature_output";                
            
            Fred_solve(Fred_equation, (std::size_t)0, this->tol, out_path);
            if (this->notifications){
                std::cout << "  Done!\n";
            }
        }

        void solve_eq_simple_iter(Class_Fred<T>& Fred_equation){

            this->check_folder(this->output_folder);
            std::string out_path;

            out_path = this->output_folder + this->file_name;

            if (this->notifications){
                std::cout << file_name << ": \t";
            }

            out_path += "_Fred_simple_iter_output";                
            
            Fred_solve(Fred_equation, (std::size_t)1, this->tol, out_path);
            if (this->notifications){
                std::cout << "  Done!\n";
            }
        }


        void solve_eq_singular(Class_Fred<T>& Fred_equation){

            this->check_folder(this->output_folder);
            std::string out_path;

            out_path = this->output_folder + this->file_name;

            if (this->notifications){
                std::cout << file_name << ": \t";
            }

            out_path += "_Fred_singular_output";                
            
            Fred_solve(Fred_equation, (std::size_t)2, this->tol, out_path);
            if (this->notifications){
                std::cout << "  Done!\n";
            }
        }


};