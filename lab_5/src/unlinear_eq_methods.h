#pragma once

#include <iostream>
#include <vector>

#define eps 1e-9

using namespace std;

// class unlinear_eq_methods
// {
// public:

//     template<typename T, typename Func>
//     vector<pair<T,T>> localization(Func foo, pair<T,T> section, size_t n);

//     template<typename T, typename Func>
//     T bisection(Func foo,  pair<T,T> section, size_t& iter_counter);
    
//     template<typename T, typename Func>
//     T chord_method(Func foo,  pair<T,T> section, size_t& iter_counter);

//     template<typename T, typename Func>
//     T chord(Func foo,  pair<T,T> section);

//     template<typename T, typename Func, typename Diff>
//     T newton_method(Func foo, Diff diff_foo, pair<T,T> section, size_t& iter_counter);

//     template<typename T, typename Func, typename J>
//     vector<T> newton_sys_method(Func foo1, Func foo2, J jacobi, const vector<T>& x0);

//     template<typename T, typename Func, typename J>
//     vector<T> newton_sys_method_mod(Func f1, Func f2, J jacobi, pair<T, T> section1, pair<T, T> section2, const vector<T>& x0);

// };


template<typename T, typename Func>
vector<pair<T,T>> localization(Func foo, pair<T,T> section, size_t n);

template<typename T, typename Func>
T bisection(Func foo,  pair<T,T> section, size_t& iter_counter);

template<typename T, typename Func>
T chord_method(Func foo,  pair<T,T> section, size_t& iter_counter);

template<typename T, typename Func>
T chord(Func foo,  pair<T,T> section);

template<typename T, typename Func, typename Diff>
T newton_method(Func foo, Diff diff_foo, pair<T,T> section, size_t& iter_counter);

template<typename T, typename Func, typename J>
vector<T> newton_sys_method(Func foo1, Func foo2, J jacobi, const vector<T>& x0);

template<typename T, typename Func, typename J>
vector<T> newton_sys_method_mod(Func f1, Func f2, J jacobi, pair<T, T> section1, pair<T, T> section2, const vector<T>& x0);