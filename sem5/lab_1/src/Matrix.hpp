#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>

#define eps 1e-16


template<class T>
class Matrix{

    private:

        

    public:

        std::size_t n;
        std::size_t rows;
        std::size_t cols;

        std::vector<std::vector<T>> array; 

        //Constructors
        Matrix<T>();
        Matrix<T>(std::size_t n);
        Matrix<T>(std::size_t rows, std::size_t cols);
        Matrix<T>(std::vector<std::vector<T>> const &array);

        //Destructor
        ~Matrix<T>(){};

        //getters
        std::size_t get_n() const;
        std::size_t get_rows() const;
        std::size_t get_cols() const;
        std::vector<std::vector<T>> get_array() const;

        //Matrix operations
        Matrix<T> add(const Matrix<T>& another_matrix) const;
        Matrix<T> subtract(const Matrix<T>& another_matrix) const;
        Matrix<T> dot(const Matrix<T>& another_matrix) const;
        Matrix<T> transpose() const;

        Matrix<T> &add_this(const Matrix<T>& another_matrix);
        Matrix<T> &subtract_this(const Matrix<T>& another_matrix);
        Matrix<T> &dot_this(const Matrix<T>& another_matrix);
        Matrix<T> &transpose_this();

        //Value operations
        Matrix<T> add(const T value) const;
        Matrix<T> subtract(const T value) const;
        Matrix<T> multiply(const T value) const;
        Matrix<T> divide(const T value) const;

        Matrix<T> &add_this(const T value);
        Matrix<T> &subtract_this(const T value);
        Matrix<T> &multiply_this(const T value);
        Matrix<T> &divide_this(const T value);



        void read_Square_Matrix(const std::string path);
        void read_Sys(const std::string path);

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


// Getters

template<class T>
std::size_t Matrix<T>::get_n() const{
    if(this->n == 0)
        std::cout << "this matrix is not squared!";

    return this->n;
}

template<class T>
std::size_t Matrix<T>::get_rows() const{
    if(this->rows == 0)
        std::cout << "null mathrix!";

    return this->rows;
}

template<class T>
std::size_t Matrix<T>::get_cols() const{
    if(this->cols == 0)
        std::cout << "null mathrix!";
    return this->cols;
}

template<class T>
std::vector<std::vector<T>> Matrix<T>::get_array() const{
    if(this-> rows || this->cols == 0)
        throw std::invalid_argument("null matrix");
    
    return this->array;
}


//Matrix operations

template<class T>
Matrix<T> Matrix<T>::add(const Matrix<T>& another_matrix) const{

    if(this->rows == 0 || this->cols == 0 || another_matrix.rows == 0 || another_matrix.cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->rows != another_matrix.rows && this->cols != another_matrix.cols)
        throw std::invalid_argument("These matrixs cannot be added first.cols != second.rows && first.rows != second.cols");
    
    Matrix result_matrix(this->array);

    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = result_matrix[i][j];
            tmp += another_matrix.array[i][j];
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            result_matrix.array[i][j] = tmp;
        }
    }
    

    return result_matrix;
}

template<class T>
Matrix<T> Matrix<T>::subtract(const Matrix<T>& another_matrix) const{

    if(this->rows == 0 || this->cols == 0 || another_matrix.rows == 0 || another_matrix.cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->rows != another_matrix.rows && this->cols != another_matrix.cols)
        throw std::invalid_argument("These matrixs cannot be substructed first.cols != second.rows && first.rows != second.cols");
    
    Matrix result_matrix(this->array);

    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = result_matrix[i][j];
            tmp -= another_matrix.array[i][j];
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            result_matrix.array[i][j] = tmp;
        }
    }

    return result_matrix;
}

template<class T>
Matrix<T> Matrix<T>::dot(const Matrix<T>& another_matrix) const{

    if(this->rows == 0 || this->cols == 0 || another_matrix.rows == 0 || another_matrix.cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->cols != another_matrix.rows){
        throw std::invalid_argument("These matrixs cannot be multiplied first.cols != second.rows");
    }

    if(this->cols == 0 || this->rows == 0 || another_matrix.rows == 0 || another_matrix.cols == 0){
        throw std::invalid_argument("These matrixs cannot be multiplied null matrix");
    }

    Matrix result_matrix(this->rows, another_matrix.cols);
    
    T tmp = 0;
    for(std::size_t i = 0; i < result_matrix.rows; ++i){
        
        for(std::size_t j = 0; j < result_matrix.cols; ++j){
        
            tmp = 0;

            for(std::size_t k = 0; k < result_matrix.cols; ++k)
                tmp += this->array[i][k] * another_matrix[k][j];

            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);

            result_matrix.array[i][j] = tmp;
        }
    } 

    return result_matrix;
}

template<class T>
Matrix<T> Matrix<T>::transpose() const{

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");


    Matrix result_matrix(this->cols, this->rows);

    for (size_t i = 0; i < n; ++i)
		for (size_t j = i + 1; j < n; ++j){
			result_matrix.array[i][j] = this->array[j][i];
			result_matrix.array[j][i] = this->array[i][j];
		}
    
    return result_matrix;
}


//Value operations
template<class T>
Matrix<T> Matrix<T>::add(const T value) const{

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix result_matrix(this->rows, this->cols);

    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j] + value;
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            result_matrix.array[i][j] = tmp;
        }
    }    

    return result_matrix;
}

template<class T>
Matrix<T> Matrix<T>::subtract(const T value) const{

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix result_matrix(this->rows, this->cols);

    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j] - value;
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            result_matrix.array[i][j] = tmp;
        }
    }    

    return result_matrix;
}


template<class T>
Matrix<T> Matrix<T>::multiply(const T value) const{

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");
    
    Matrix result_matrix(this->rows, this->cols);

    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j] * value;
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            result_matrix.array[i][j] = tmp;
        }
    }    

    return result_matrix;
}


template<class T>
Matrix<T> Matrix<T>::divide(const T value) const{

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix result_matrix(this->rows, this->cols);

    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j] / value;
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            result_matrix.array[i][j] = tmp;
        }
    }    

    return result_matrix;
}









template<class T>
Matrix<T> &Matrix<T>::add_this(const Matrix<T>& another_matrix){

    if(this->rows == 0 || this->cols == 0 || another_matrix.rows == 0 || another_matrix.cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->rows != another_matrix.rows && this->cols != another_matrix.cols)
        throw std::invalid_argument("These matrixs cannot be added first.cols != second.rows && first.rows != second.cols");
    
    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j];
            tmp += another_matrix.array[i][j];
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            this->array[i][j] = tmp;
        }
    }
    

    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::subtract_this(const Matrix<T>& another_matrix){

    if(this->rows == 0 || this->cols == 0 || another_matrix.rows == 0 || another_matrix.cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->rows != another_matrix.rows && this->cols != another_matrix.cols)
        throw std::invalid_argument("These matrixs cannot be substructed first.cols != second.rows && first.rows != second.cols");
    
    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j];
            tmp -= another_matrix.array[i][j];
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            this->array.array[i][j] = tmp;
        }
    }

    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::dot_this(const Matrix<T>& another_matrix){

    if(this->rows == 0 || this->cols == 0 || another_matrix.rows == 0 || another_matrix.cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->cols != another_matrix.rows){
        throw std::invalid_argument("These matrixs cannot be multiplied first.cols != second.rows");
    }

    if(this->cols == 0 || this->rows == 0 || another_matrix.rows == 0 || another_matrix.cols == 0){
        throw std::invalid_argument("These matrixs cannot be multiplied null matrix");
    }

    std::vector<std::vector<T>> result_array(this->rows, std::vector(another_matrix.cols));
    
    T tmp = 0;
    for(std::size_t i = 0; i < this->rows; ++i){
        
        for(std::size_t j = 0; j < another_matrix.cols; ++j){
        
            tmp = 0;

            for(std::size_t k = 0; k < another_matrix.cols; ++k)
                tmp += this->array[i][k] * another_matrix[k][j];

            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);

            result_array[i][j] = tmp;
        }
        
        this->cols = another_matrix.cols;
        this->array = result_array;
    } 

    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::transpose_this(){

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");


    if(this->rows != this->cols){
        Matrix result_matrix(this->cols, this->rows);

        std::vector<std::vector<T>> result_array(this->cols, std::vector<T>(this->rows));

        for (std::size_t i = 0; i < this->rows; ++i)
            for (std::size_t j = i; j < this->cols; ++j){
                result_array[j][i] = this->array[i][j];
                result_array[i][j] = this->array[j][i];
            }
        this->rows = result_array.size();
        this->cols = result_array[0].size();
        this->array = result_array;

    }
    else{
        T tmp;
        for (std::size_t i = 0; i < this->rows; ++i)
            for (std::size_t j = i+1; j < this->cols; ++j){
                tmp = this->array[i][j];
                this->array[i][j] = this->array[j][i];
                this->array[j][i] = tmp;
            }
    }    

    return *this;
}


//Value operations
template<class T>
Matrix<T> &Matrix<T>::add_this(const T value){

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;

    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j] + value;        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            this->array[i][j] = tmp;
        }
    }    

    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::subtract_this(const T value){

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j] - value;
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            this->array[i][j] = tmp;
        }
    }    

    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::multiply_this(const T value){

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");
    
    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j] * value;
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            this->array[i][j] = tmp;
        }
    }    

    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::divide_this(const T value){

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;
    for(std::size_t i = 0; i < this->n; ++i){        
        for(std::size_t j = 0; j < this->n; ++j){

            tmp = this->array[i][j] / value;
        
            if(fabs(tmp) < eps)
                tmp *= (1 / tmp > 0) ? (1) : (-1);
            
            this->array[i][j] = tmp;
        }
    }    

    return *this;
}












template<class T>
void Matrix<T>::read_Square_Matrix(const std::string path){

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

    if(this->rows == 0 || this->cols == 0)
        throw std::invalid_argument("null matrix");

    for (std::size_t i = 0; i < this->rows; ++i) {
        for (std::size_t j = 0; j < this->cols; ++j) 
            std::cout << setw(15) << this->array[i][j];		
        std::cout << "\n";
    }
    std::cout << "\n";
}