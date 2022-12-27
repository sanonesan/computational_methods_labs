#pragma once
#include <cmath>
#include <vector>

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;


template<typename T>
T test0(T x);

template<typename T>
T test0_diff(T x);

template<typename T>
T test1(T x);

template<typename T>
T test1_diff(T x);

template<typename T>
T test3(T x);

template<typename T>
T test3_diff(T x);

template<typename T>
T test4(T x);

template<typename T>
T test4_diff(T x);

template<typename T>
T test4f1(T x, T y);

template<typename T>
T test4f2(T x, T y);

template<typename T>
matrix<T> test5Jacobi(T x, T y);


template<typename T>
T test5f1(T x, T y);

template<typename T>
T test5f2(T x, T y);

template<typename T>
matrix<T> test5Jacobi(T x, T y);

template<typename T>
T test6(T x);

template<typename T>
T test6_diff(T x);

