#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>

#define eps 1e-16


template<class T>
class Matrix;

template<class T>
class Vector;



//----------------------------------------//
//----------------VECTOR------------------//
//----------------------------------------//


template<class T>
class Vector{


    private:
        std::size_t _size;
    public:

        friend Matrix<T>;

        std::vector<T> _vector;

        Vector<T>();
        Vector<T>(const std::size_t size);
        Vector<T>(const std::vector<T>& vec);
        Vector<T>(const Vector<T>& vec);
        ~Vector<T>(){};

        std::size_t size() const;

        Vector<T> dot(const Matrix<T>& matrix);

        //vector norms

        T norm_1() const;
        T norm_inf() const;
        T norm_euclid() const;


};
template<class T>
std::ostream& operator<<(std::ostream &out, const Vector<T>& vec);


//----------------------------------------//
//----------------MATRIX------------------//
//----------------------------------------//


template<class T>
class Matrix{


    private:

        std::size_t _rows;
        std::size_t _cols;


    public:

        friend Vector<T>;

        std::vector<std::vector<T>> _array; 
        
        //Constructors
        Matrix<T>();
        Matrix<T>(std::size_t _n);
        Matrix<T>(std::size_t _rows, std::size_t _cols);
        Matrix<T>(std::vector<std::vector<T>> const &_array);
        //Copy constructor
        Matrix<T>(const Matrix<T>& matrix);
        //Destructor
        ~Matrix<T>(){};

        //getters
        std::size_t get_n() const;
        std::size_t get_rows() const;
        std::size_t get_cols() const;
        std::vector<std::vector<T>> get_array() const;

        //-----------RESULT = NEW_MATRIX-----------//

        //Matrix operations


        Matrix<T> add(const Matrix<T>& another_matrix) const;
        Matrix<T> subtract(const Matrix<T>& another_matrix) const;
        Matrix<T> dot(const Matrix<T>& another_matrix) const;
        Matrix<T> transpose() const;


        //Value operations

        Matrix<T> add(const T value) const;
        Matrix<T> subtract(const T value) const;
        Matrix<T> multiply(const T value) const;
        Matrix<T> divide(const T value) const;


        //--------INFLUENCE_CURRENT_MATRIX---------//

        //Matrix operations

        Matrix<T> &add_this(const Matrix<T>& another_matrix);
        Matrix<T> &subtract_this(const Matrix<T>& another_matrix);
        Matrix<T> &dot_this(const Matrix<T>& another_matrix);
        Matrix<T> &transpose_this();


        //Value operations

        Matrix<T> &add_this(const T value);
        Matrix<T> &subtract_this(const T value);
        Matrix<T> &multiply_this(const T value);
        Matrix<T> &divide_this(const T value);


        //matrix norms

        T norm_1() const;
        T norm_inf() const;

        void read_Square_Matrix(const std::string path);
        void read_Matrix(const std::string path);
        void print();


        //Bool operations

        bool operator==(const Matrix<T>& another_matrix);
        bool operator!=(const Matrix<T>& manother_matrix);


        


        //Matrix operations

        Matrix<T> operator+(const Matrix<T>& another_matrix);
        Matrix<T> operator-(const Matrix<T>& another_matrix);
        Matrix<T> operator*(const Matrix<T>& another_matrix);


        //Value operations

        Matrix<T> operator+(const T value);
        Matrix<T> operator-(const T value);
        Matrix<T> operator*(const T value);

        //Matrix operations

        Matrix<T> &operator=(const Matrix<T>& another_matrix);
        Matrix<T> operator+=(const Matrix<T>& another_matrix);
        Matrix<T> operator-=(const Matrix<T>& another_matrix);
        Matrix<T> operator*=(const Matrix<T>& another_matrix);

        //Value operations

        Matrix<T> operator+=(const T value);
        Matrix<T> operator-=(const T value);
        Matrix<T> operator*=(const T value);

        Vector<T> dot(const Vector<T>& vector);
        Vector<T> dot(const std::vector<T>& vector);
        //friend Matrix<T> dot_thins(const Vector<T>& vector);

};

template<class T>
std::ostream& operator<<(std::ostream &out, const Matrix<T>& matrix);



//----------------------------------------//
//----------------MATRIX------------------//
//----------------------------------------//


template<class T>
Matrix<T>::Matrix(){
    this->_rows = 0;
    this->_cols = 0;
}


template<class T>
Matrix<T>::Matrix(std::size_t _n){
    this->_cols = _n;
    this->_rows = _n;
    this->_array = std::vector<std::vector<T>> (_n, std::vector<T>(_n));
}


template<class T>
Matrix<T>::Matrix(std::size_t _rows, std::size_t _cols){

    this->_rows = _rows;
    this->_cols = _cols;
    this->_array = std::vector<std::vector<T>> (_rows, std::vector<T>(_cols));

    for(std::size_t i = 0; i < this->_rows; ++i){
        for(std::size_t j = 0; j < this->_cols; ++j)
            this->_array[i][j] = 0.;
    }
}


template<class T>
Matrix<T>::Matrix(std::vector<std::vector<T>> const &_array){

    if(_array.size() == 0 || _array[0].size() == 0)
        throw std::invalid_argument("Size of _array should be > 0");

    this->_rows = _array.size();
    this->_cols = _array[0].size();
    this->_array = _array;
}

//Copy constructor
template<class T>
Matrix<T>::Matrix(const Matrix<T>& matrix){
    this->_rows = matrix._rows;
    this->_cols = matrix._cols;
    this->_array = matrix._array;
}


// Getters

template<class T>
std::size_t Matrix<T>::get_rows() const{
    if(this->_rows == 0)
        std::cout << "null matrix!";

    return this->_rows;
}

template<class T>
std::size_t Matrix<T>::get_cols() const{
    if(this->_cols == 0)
        std::cout << "null matrix!";
    return this->_cols;
}

template<class T>
std::vector<std::vector<T>> Matrix<T>::get_array() const{
    if(this-> _rows || this->_cols == 0)
        throw std::invalid_argument("null matrix");
    
    return this->_array;
}



/*
/////////////////////////////////////////////
//-----------RESULT = NEW_MATRIX-----------//
/////////////////////////////////////////////
*/

//Matrix operations


template<class T>
Matrix<T> Matrix<T>::add(const Matrix<T>& another_matrix) const{

    if(this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->_rows != another_matrix._rows && this->_cols != another_matrix._cols)
        throw std::invalid_argument("These matrixs cannot be added first._cols != second._rows && first._rows != second._cols");
    
    Matrix<T> result_matrix(this->_array);

    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = result_matrix[i][j];
            tmp += another_matrix._array[i][j];
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            result_matrix._array[i][j] = tmp;
        }
    }
    

    return result_matrix;
}

template<class T>
Matrix<T> Matrix<T>::subtract(const Matrix<T>& another_matrix) const{

    if(this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->_rows != another_matrix._rows && this->_cols != another_matrix._cols)
        throw std::invalid_argument("These matrixs cannot be substructed first._cols != second._rows && first._rows != second._cols");
    
    Matrix<T> result_matrix(this->_array);

    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = result_matrix[i][j];
            tmp -= another_matrix._array[i][j];
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            result_matrix._array[i][j] = tmp;
        }
    }

    return result_matrix;
}

template<class T>
Matrix<T> Matrix<T>::dot(const Matrix<T>& another_matrix) const{

    if(this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->_cols != another_matrix._rows){
        throw std::invalid_argument("These matrixs cannot be multiplied first._cols != second._rows");
    }

    if(this->_cols == 0 || this->_rows == 0 || another_matrix._rows == 0 || another_matrix._cols == 0){
        throw std::invalid_argument("These matrixs cannot be multiplied null matrix");
    }

    Matrix<T> result_matrix(this->_rows, another_matrix._cols);
    
    T tmp = 0;
    for(std::size_t i = 0; i < result_matrix._rows; ++i){
        
        for(std::size_t j = 0; j < result_matrix._cols; ++j){
        
            tmp = 0;

            for(std::size_t k = 0; k < result_matrix._cols; ++k)
                tmp += this->_array[i][k] * another_matrix._array[k][j];

            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_matrix._array[i][j] = tmp;
        }
    } 

    return result_matrix;
}

template<class T>
Matrix<T> Matrix<T>::transpose() const{

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");


    Matrix<T> result_matrix(this->_cols, this->_rows);

    for (std::size_t i = 0; i < this->_rows; ++i)
        for (std::size_t j = i; j < this->_cols; ++j){
            result_matrix._array[i][j] = this->_array[j][i];
            result_matrix._array[j][i] = this->_array[i][j];
        }
    
    return result_matrix;
}


//Value operations
template<class T>
Matrix<T> Matrix<T>::add(const T value) const{

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix<T> result_matrix(this->_rows, this->_cols);

    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j] + value;
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            result_matrix._array[i][j] = tmp;
        }
    }    

    return result_matrix;
}

template<class T>
Matrix<T> Matrix<T>::subtract(const T value) const{

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix<T> result_matrix(this->_rows, this->_cols);

    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j] - value;
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            result_matrix._array[i][j] = tmp;
        }
    }    

    return result_matrix;
}


template<class T>
Matrix<T> Matrix<T>::multiply(const T value) const{

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");
    
    Matrix<T> result_matrix(this->_rows, this->_cols);

    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j] * value;
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            result_matrix._array[i][j] = tmp;
        }
    }    

    return result_matrix;
}


template<class T>
Matrix<T> Matrix<T>::divide(const T value) const{

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    Matrix<T> result_matrix(this->_rows, this->_cols);

    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j] / value;
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            result_matrix._array[i][j] = tmp;
        }
    }    

    return result_matrix;
}


/*
/////////////////////////////////////////////
//--------INFLUENCE_CURRENT_MATRIX---------//
/////////////////////////////////////////////
*/

//Matrix operations

template<class T>
Matrix<T> &Matrix<T>::add_this(const Matrix<T>& another_matrix){

    if(this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->_rows != another_matrix._rows && this->_cols != another_matrix._cols)
        throw std::invalid_argument("These matrixs cannot be added first._cols != second._rows && first._rows != second._cols");
    
    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j];
            tmp += another_matrix._array[i][j];
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            this->_array[i][j] = tmp;
        }
    }
    

    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::subtract_this(const Matrix<T>& another_matrix){

    if(this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->_rows != another_matrix._rows && this->_cols != another_matrix._cols)
        throw std::invalid_argument("These matrixs cannot be substructed first._cols != second._rows && first._rows != second._cols");
    
    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j];
            tmp -= another_matrix._array[i][j];
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            this->_array[i][j] = tmp;
        }
    }

    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::dot_this(const Matrix<T>& another_matrix){

    if(this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->_cols != another_matrix._rows){
        throw std::invalid_argument("These matrixs cannot be multiplied first._cols != second._rows");
    }

    if(this->_cols == 0 || this->_rows == 0 || another_matrix._rows == 0 || another_matrix._cols == 0){
        throw std::invalid_argument("These matrixs cannot be multiplied null matrix");
    }

    std::vector<std::vector<T>> result_array(this->_rows, std::vector<T> (another_matrix._cols));

    T tmp = 0;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < another_matrix._cols; ++j){
        
            tmp = 0.;

            for(std::size_t k = 0; k < another_matrix._cols; ++k)
                tmp += this->_array[i][k] * another_matrix._array[k][j];

            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

            result_array[i][j] = tmp;
        }
    } 

    this->_cols = another_matrix._cols;
    this->_array = result_array;
    
    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::transpose_this(){

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");


    if(this->_rows != this->_cols){
        Matrix result_matrix(this->_cols, this->_rows);

        std::vector<std::vector<T>> result_array(this->_cols, std::vector<T>(this->_rows));

        for (std::size_t i = 0; i < this->_rows; ++i)
            for (std::size_t j = i; j < this->_cols; ++j){
                result_array[j][i] = this->_array[i][j];
                result_array[i][j] = this->_array[j][i];
            }
        this->_rows = result_array.size();
        this->_cols = result_array[0].size();
        this->_array = result_array;

    }
    else{
        T tmp;
        for (std::size_t i = 0; i < this->_rows; ++i)
            for (std::size_t j = i+1; j < this->_cols; ++j){
                tmp = this->_array[i][j];
                this->_array[i][j] = this->_array[j][i];
                this->_array[j][i] = tmp;
            }
    }    

    return *this;
}


//Value operations


template<class T>
Matrix<T> &Matrix<T>::add_this(const T value){

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;

    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j] + value;        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            this->_array[i][j] = tmp;
        }
    }    

    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::subtract_this(const T value){

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j] - value;
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            this->_array[i][j] = tmp;
        }
    }    

    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::multiply_this(const T value){

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");
    
    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j] * value;
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            this->_array[i][j] = tmp;
        }
    }    

    return *this;
}


template<class T>
Matrix<T> &Matrix<T>::divide_this(const T value){

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    T tmp;
    for(std::size_t i = 0; i < this->_rows; ++i){        
        for(std::size_t j = 0; j < this->_cols; ++j){

            tmp = this->_array[i][j] / value;
        
            if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;
            
            this->_array[i][j] = tmp;
        }
    }    

    return *this;
}



template<class T>
T Matrix<T>::norm_1() const{
    T norm = 0.;
    T tmp = 0.;
    for(std::size_t i = 0; i < this->_rows; ++i){
        tmp = 0.;
        for(std::size_t j = 0; j < this->_cols; ++j){
            tmp += fabs(this->_array[j][i]);
        }
        if (norm < tmp) norm = tmp;
    }
    if(fabs(norm) < eps && norm != 0.){
        norm = 0.;
        if (1 / norm < 0) norm *= -1;
    }
    return norm;
}


template<class T>
T Matrix<T>::norm_inf() const{
    T norm = 0.;
    T tmp = 0.;
    for(std::size_t i = 0; i < this->_rows; ++i){
        tmp = 0.;
        for(std::size_t j = 0; j < this->_cols; ++j){
            tmp += fabs(this->_array[i][j]);
        }
        if (norm < tmp) norm = tmp;
    }
    if(fabs(norm) < eps && norm != 0.){
        norm = 0.;
        if (1 / norm < 0) norm *= -1;
    }
    return norm;
}


template<class T>
void Matrix<T>::read_Square_Matrix(const std::string path){

    std::ifstream fin;
    fin.open(path);
    if(!fin.is_open()){
        throw std::invalid_argument("Smth's wrong with path");
    }
    fin >> this->_rows;
    this->_cols = this->_rows;
    fin.get();    
    std::vector<std::vector<T>> _array(this->_rows, std::vector<T>(this->_cols));

    for(std::size_t i = 0; i < this->_rows; ++i){
        for(std::size_t j = 0; j < this->_cols; ++j){
            fin >> _array[i][j];
        }
    }

    fin.close();
    this->_array = _array;
}


template<class T>
void Matrix<T>::read_Matrix(const std::string path){

    std::ifstream fin;
    fin.open(path);
    if(!fin.is_open()){
        throw std::invalid_argument("Smth's wrong with path");
    }
    fin >> this->_rows;
    fin >> this->_cols;
    fin.get();    
    std::vector<std::vector<T>> colsarray(this->_rows, std::vector<T>(this->_cols));

    for(std::size_t i = 0; i < this->_rows; ++i){
        for(std::size_t j = 0; j < this->_cols; ++j){
            fin >> _array[i][j];
        }
    }

    fin.close();
    this->_array = _array;
}


template<class T>
void Matrix<T>::print(){

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    for (std::size_t i = 0; i < this->_rows; ++i) {
        for (std::size_t j = 0; j < this->_cols; ++j) 
            std::cout << std::setw(15) << this->_array[i][j];		
        std::cout << "\n";
    }
    std::cout << "\n";
}


template<class T>
bool Matrix<T>::operator==(const Matrix<T>& another_matrix){
    if(this->_rows == 0 || this->_cols == 0 || another_matrix._rows == 0 || another_matrix._cols == 0)
        return false;

    if(this->_rows != another_matrix._rows || this->_cols != another_matrix._cols)
        return false;

    for(std::size_t i = 0; i < this->_rows; ++i){
        for(std::size_t j = 0; j < this -> _cols; ++j){
            if (fabs(this->_array[i][j] - another_matrix._array[i][j]) > eps) 
                return false;
        }
    }

    return true;
    
}

template<class T>
bool Matrix<T>::operator!=(const Matrix<T>& another_matrix){
    return !operator==(another_matrix);
}


//Matrix operations;


template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& another_matrix){    
    return this->add(another_matrix);
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& another_matrix){
    return this->subtract(another_matrix);
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& another_matrix){
    return this->dot(another_matrix);
}


//Value operations

template<class T>
Matrix<T> Matrix<T>::operator+(const T value){    
    return this->add(value);
}

template<class T>
Matrix<T> Matrix<T>::operator-(const T value){
    return this->subtract(value);
}

template<class T>
Matrix<T> Matrix<T>::operator*(const T value){
    return this->multiply(value);
}


//Matrix operations;
template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T>& another_matrix){   
    this->_rows = another_matrix._rows;
    this->_cols = another_matrix._cols;
    this->_array = another_matrix._array;
    return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator+=(const Matrix<T>& another_matrix){    
    return this->add_this(another_matrix);
}

template<class T>
Matrix<T> Matrix<T>::operator-=(const Matrix<T>& another_matrix){
    return this->subtract_this(another_matrix);
}

template<class T>
Matrix<T> Matrix<T>::operator*=(const Matrix<T>& another_matrix){
    return this->dot_this(another_matrix);
}


//Value operations

template<class T>
Matrix<T> Matrix<T>::operator+=(const T value){    
    return this->add_this(value);
}

template<class T>
Matrix<T> Matrix<T>::operator-=(const T value){
    return this->subtract_this(value);
}

template<class T>
Matrix<T> Matrix<T>::operator*=(const T value){
    return this->multiply_this(value);
}


//iostream 

template<class T>
std::ostream& operator<<(std::ostream &out, const Matrix<T>& matrix){
    if(matrix.get_rows() == 0 || matrix.get_cols() == 0)
        throw std::invalid_argument("null matrix");

    for (std::size_t i = 0; i < matrix.get_rows(); ++i) {
        for (std::size_t j = 0; j < matrix.get_cols(); ++j) 
            out << std::setw(15) << matrix._array[i][j];		
        out << "\n";
    }
    out << "\n";

    return out;
}

template<class T>
Vector<T> Matrix<T>::dot(const Vector<T>& vector){

    if(this->_rows == 0 || this->_cols == 0)
        throw std::invalid_argument("null matrix");

    if(vector.size() == 0)
        throw std::invalid_argument("null vector");

    if(this->_cols != vector._size)
        throw std::invalid_argument("matrix.cols != vector.size");

    Vector<T> result_vector(this->_rows);

    T tmp;

    for(std::size_t i = 0; i < this->_rows; ++i){
        tmp = 0.;
        for(std::size_t j = 0; j < this->_cols; ++j){
            tmp += this->_array[i][j] * vector._vector[j];
        }
        if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

        result_vector._vector[i] = tmp;
    }

    return result_vector;
}


template<class T>
Vector<T> Matrix<T>::dot(const std::vector<T>& vector){
    Vector<T> tmp(vector);
    return this->dot(tmp);
};






//----------------------------------------//
//----------------VECTOR------------------//
//----------------------------------------//



template<class T>
Vector<T>::Vector(){
    this->_size = 0;
};


template<class T>
Vector<T>::Vector(const std::size_t size){
    this->_size = size;
    this->_vector = std::vector<T> (size);
};


template<class T>
Vector<T>::Vector(const std::vector<T>& vec){
    this->_vector = vec;
    this->_size = this->_vector.size();
};

template<class T>
Vector<T>::Vector(const Vector<T>& vec){
    this->_size = vec._size;
    this->_vector = vec._vector;
};

template<class T>
std::size_t Vector<T>::size() const{
    return this->_size;
};


template<class T>
Vector<T> Vector<T>::dot(const Matrix<T>& matrix){

    if(matrix._rows == 0 || matrix._cols == 0)
        throw std::invalid_argument("null matrix");

    if(this->_size == 0)
        throw std::invalid_argument("null vector");

    if(this->_size != matrix._rows)
        throw std::invalid_argument("vector._size != matrix._rows");

    Vector<T> result_vector(this->_size);

    T tmp;

    for(std::size_t i = 0; i < this->_size; ++i){
        tmp = 0.;
        for(std::size_t j = 0; j < matrix._rows; ++j){
            tmp += this->_vector[i] * matrix._array[i][j];
        }
        if(fabs(tmp) < eps)
                if (1 / tmp < 0) tmp *= -1;

        result_vector._vector[i] = tmp;
    }

    return result_vector;
}


/*
	Октаэдрическая норма ( || * ||_{1})
*/
template<class T>
T Vector<T>::norm_1() const{
	
	T norm = 0;
	for(std::size_t i = 0; i < this->_size; ++i)
		norm += fabs(this->_vector[i]);

	return norm < eps ? 0. : norm;
}


/*
	Кубическая норма ( || * ||_{inf})
*/
template<class T> 
T Vector<T>::norm_inf() const{
	
	T norm = fabs(this->_vector[0]);
	for(std::size_t i = 1; i < this->_size; ++i)
		if(fabs(norm - fabs(this->_vector[i]) < eps))
			norm = fabs(this->_vector[i]);
	
	return norm < eps ? 0. : norm;
}


/*
	Евклидова норма ( || * ||_{eucl})
*/
template<class T> 
T Vector<T>::norm_euclid() const{

	T norm = 0;
	for(std::size_t i = 0; i < this->_size; ++i)
		norm += fabs(vehis->_vector[i] * his->_vector[i]);

	return norm < eps ? 0. : sqrt(norm);
}


template<class T>
std::ostream& operator<<(std::ostream &out, const Vector<T>& vec){

    out << "( ";
    for(std::size_t i = 0; i < vec.size() - 1; ++i){
        out << vec._vector[i] << "\t";
    }
    out << vec._vector[vec.size() - 1] << " )^T";
    
    return out;
}