#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>

template<class T>
class Matrix{

    private:

        std::size_t n;
        std::size_t rows;
        std::size_t cols;
        std::vector<std::vector<T>> array;        

    public:

        //default constructor
        Matrix<T>();
        Matrix<T>(std::size_t n);
        Matrix<T>(std::size_t rows, std::size_t cols);
        Matrix<T>(std::vector<std::vector<T>> const &array);


        void read_Matrix(const std::string path);
        void read_Sys(const std::string path);


        //default destructor
        ~Matrix(){};

        void Print();
        
};


template<class T>
Matrix<T>::Matrix(){
    this->n = 0;
    this->rows = 0;
    this->cols = 0;
}

template<class T>
Matrix<T>::Matrix(std::size_t _n){
    this->n = _n;
    this->cols = _n;
    this->rows = _n;
    this->array = std::vector<std::vector<T>> (_n, std::vector<T>(_n));
}

template<class T>
Matrix<T>::Matrix(std::size_t _rows, std::size_t _cols){

    if(_rows == _cols){ this->n = _rows; }
    else{ this->n = 0; }
        
    this->rows = _rows;
    this->cols = _cols;
    this->array = std::vector<std::vector<T>> (_rows, std::vector<T>(_cols));

    for(std::size_t i = 0; i < this->n; ++i){
        for(std::size_t j = 0; j < this->n; ++j)
            this->array[i][j] = 0.;
    }
}

template<class T>
Matrix<T>::Matrix(std::vector<std::vector<T>> const &array){

    if(array.size() == 0 || array[0].size() == 0)
        throw std::invalid_argument("Size of array should be > 0");

    this->rows = array.size();
    this->cols = array[0].size();
    this->array = array;
}

template<class T>
void Matrix<T>::read_Matrix(const std::string path){

    std::ifstream fin;
    fin.open(path);
    if(!fin.is_open()){
        throw std::invalid_argument("Smth's wrong with path");
    }

    fin >> this->n;
    fin.get();    
    std::vector<std::vector<T>> array(this->n, std::vector<T>(this->n));

    for(std::size_t i = 0; i < this->n; ++i){
        for(std::size_t j = 0; j < this->n; ++j){
            fin >> array[i][j];
        }
    }

    fin.close();

    this->rows = this->n;
    this->cols = this->n;
    this->array = array;
}



template<typename T>
void Matrix<T>::Print(){
    for (std::size_t i = 0; i < this->n; ++i) {
        for (std::size_t j = 0; j < this->n; ++j) 
            std::cout << setw(15) << this->array[i][j];		
        std::cout << "\n";
    }
    std::cout << "\n";
}