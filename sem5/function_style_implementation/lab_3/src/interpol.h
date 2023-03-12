#pragma once

#include <iostream>
#include <fstream>

using namespace std;

#define eps 1e-16

double interpol_lagrange(double x, double*& x_mas, double*& y_mas, size_t& n);