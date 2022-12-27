#include <cmath>
#include <vector>
#include "unlinear_methods_test_functions.h"

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;


template<typename T>
T test0(T x)
{
	return (x - 1) * (x - 1) * (x - 3) * (x - 3);
}

template<typename T>
T test0_diff(T x)
{
	return 4 * (-6 + 11 * x - 6 * x * x + x * x * x);
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
pair<T, T> test5(T x, T y)
{
	return make_pair(x * x + y * y + x + y - 8, x * x + y * y + x * y - 7);
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
