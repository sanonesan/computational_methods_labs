#include <cmath>
#include <vector>
#include "unlinear_methods_test_functions.h"

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;

template<typename T>
T test0(T x)
{
	return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

template<typename T>
T test0_diff(T x)
{
	return 0.121495 - 1.5119 * x + 5.9535 * x * x - 9.28 * x * x * x + 5. * x * x * x * x;
}


template<typename T>
T test1(T x)
{
	return sqrt(x + 1) - 1;
}

template<typename T>
T test1_diff(T x)
{
	return 1. / (2. * sqrt(1. + x));
}

template<typename T>
T test3(T x)
{
	return 3 * x * x * x * x - 5 * x * x * x + x;
}

template<typename T>
T test3_diff(T x)
{
	return 3 * 4 * x * x * x - 5 * 3 * x * x + 1;
}

template<typename T>
T test4(T x)
{
	return 35 * x * x * x - 67 * x * x - 3 * x + 3;
}

template<typename T>
T test4_diff(T x)
{
	return -3. - 134. * x + 105. * x * x;
}





template<typename T>
T test4f1(T x, T y)
{
	return x * x - y * y - 15.;
}

template<typename T>
T test4f2(T x, T y)
{
	return x * y + 4.;
}

template<typename T>
matrix<T> test4Jacobi(T x, T y)
{
	return matrix<T>{
		{2. * x, -2. * y},
		{ y, x }
	};
}








template<typename T>
T test5f1(T x, T y)
{
	return x * x + y * y + x + y - 8.;
}

template<typename T>
T test5f2(T x, T y)
{
	return x * x + y * y + x * y - 7.;
}

template<typename T>
matrix<T> test5Jacobi(T x, T y)
{
	return matrix<T>{
		{1. + 2. * x, 1. + 2. * y},
		{ 2. * x + y, x + 2. * y }
	};
}


template<typename T>
T test6(T x) {
	return 7 * sin(x) - 4;
}
template<typename T>
T test6_diff(T x) {
	return 7 * cos(x);
}


template<typename T>
T test7(T x){
	return pow(x-1, 7);
}

template<typename T>
T test7_diff(T x){
	return 7 * pow(x-1, 6);
}