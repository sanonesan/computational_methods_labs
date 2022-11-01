#pragma once
#define eps 1e-4
#include <string>
//#include <cstddef>
#include <iostream>


using namespace std;

double** readMatrixFromFile(string path, size_t& n);

void print_Sys_Matrix(double**& A, size_t& n);

void print_Matrix(double**& A, size_t& n);

void print_Vector(double* &vec, size_t& n);

int copy_Matrix_A_to_B(double**& A, double**& B, size_t& n);

int copy_Sys_A_to_B(double**& A, double**& B, size_t& n);

int copy_Vector_A_to_B(double*& A, double*& B, size_t& n);

int matrix_make_E(double** &A, size_t &n);

int matrix_number_mult(double**& A, double& w, size_t &n);

int matrix_mult(double** &A, double** &B, double**& T, size_t n);

int matrix_vector_mult(double** A, double* B, double*& T, size_t n);

int normalize_Sys_MatrixLine(size_t& k, double**& A, size_t& n);

int check_matrix_zero(double**& A, size_t& n);

int check_sys_matrix_zero(double**& A, size_t& n);

int check_vector_zero(double*& A, size_t& n);



//-------------//

int norm_1_vector(double*& A, double& x, size_t& n);

int norm_inf_vector(double*& A, double& x, size_t& n);

int norm_1_matrix(double**& A, double& x, size_t& n);

int norm_inf_matrix(double**& A, double& x, size_t& n);

int discrepancy(double**& A, double*& x, double& d, size_t& n, int flag);

//int discrepancy(double**& A, double*& x, double& d1, double& dinf, size_t& n);

int transpose_matrix(double**& A, size_t& n);

int reverse_matrix(double**& A, double**& B, size_t& n);

int cond_1_matrix(double**& A, double& cond1, size_t& n);

int cond_inf_matrix(double**& A, double& condinf, size_t& n);

int vector_A_min_B(double*& A, double* B, double*& res, size_t& n);

int vector_Sys_A_min_B(double**& A, double** B, double*& res, size_t& n);

int vector_valuation(double*& X, double*& _X, double& res, size_t& n, int flag);

int vector_valuation_Sys(double**& X, double**& _X, double& res, size_t& n, int flag);

int vozm(double**& A, double*& x_save, double& max1, double& max_inf, size_t& n);

int cond_matrix(double**& A, double& cond1, size_t& n, int flag);

int reverse_matrix_qr(double**& A, double**& B, size_t& n);

int LDU_method(double**& A, double**& L, double**& D, double**& U, size_t& n);