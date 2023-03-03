#pragma once

#include<vector>
#include<cmath>

#define eps 1e-16

using namespace std;


/* 
    template<typename T>
    T norm_1(const vector<T> &vec);

    template<typename T>
    T norm_inf(const vector<T> &vec);

    template<typename T>
    T norm_euclid(const vector<T> &vec);
*/

template<typename T>
void print_vec(vector<T> vec){
	cout << "\n( ";
	for(size_t i = 0; i < vec.size() - 1; ++i)
		cout <<  vec[i] << "\t";
	cout << vec[vec.size() - 1] << " )^T \n";
}

/*
	Октаэдрическая норма ( || * ||_{1})
*/
template<typename T>
T norm_1(const vector<T> &vec){
	
	T res = 0;
	for(size_t i = 0; i < vec.size(); ++i)
		res += fabs(vec[i]);

	return res < eps ? 0. : res;
}


/*
	Кубическая норма ( || * ||_{inf})
*/
template<typename T>
T norm_inf(const vector<T> &vec){
	
	T res = fabs(vec[0]);
	for(size_t i = 1; i < vec.size(); ++i)
		if(fabs(res - fabs(vec[i]) < eps))
			res = fabs(vec[i]);
	
	return res < eps ? 0. : res;
}


/*
	Евклидова норма ( || * ||_{eucl})
*/
template<typename T>
T norm_euclid(const vector<T> &vec){

	T res = 0;
	for(size_t i = 0; i < vec.size(); ++i)
		res += fabs(vec[i] * vec[i]);

	return res < eps ? 0. : sqrt(res);
}