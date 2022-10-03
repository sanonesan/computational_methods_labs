#pragma once

#include <cstddef>

using namespace std;
//--------------------//
int method_gauss(double**& A, size_t& n, double*& x);

//subfunctions1
int straight_course(double**& A, size_t& n);
int reverse_course(double**& A, size_t& n, double*& x);

//subfunctions2
int major_element(double**& A, size_t& n);
//int normalize_Sys_MatrixLine(size_t& k, double**& A, size_t& n);
int remove_MatrixColumnElements_UnderLine(size_t& k, double**& A, size_t& n);

