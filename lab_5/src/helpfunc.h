#pragma once
#include <vector>

using namespace std;

template<typename T>
using matrix = vector<vector<T>>;



template<typename T>
void print_vec(const vector<T>& vec);

template<typename T>
void print_vec(const vector<pair<T, T>>& vec);

template<typename T>
vector<T> mult(const matrix<T>& A1, const vector<T>& A2);

template<typename T>
T norm(const vector<T>& x);

template<typename T>
matrix<T> inv(const matrix<T>& m);

template<typename T>
void error_order(vector<T> &x_mas, T x);
