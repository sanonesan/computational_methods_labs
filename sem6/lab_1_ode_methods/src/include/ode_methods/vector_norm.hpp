#ifndef ODE_VECTOR_NORMS_FOR_METHODS_HPP
#define ODE_VECTOR_NORMS_FOR_METHODS_HPP

#include<iostream>
#include<vector>
#include<cmath>

#define eps 1e-16


/* 
    template<typename T>
    T norm_1(const std::vector<T> &vec);

    template<typename T>
    T norm_inf(const std::vector<T> &vec);

    template<typename T>
    T norm_euclid(const std::vector<T> &vec);
*/


template<typename T>
void print_vec(std::vector<T> vec){
	std::cout << "\n( ";
	for(std::size_t i = 0; i < vec.size() - 1; ++i)
		std::cout <<  vec[i] << "\t";
	std::cout << vec[vec.size() - 1] << " )^T \n";
}

/*
	Октаэдрическая норма ( || * ||_{1})
*/
template<typename T>
T norm_1(const std::vector<T> &vec){
	
	T res = 0;
	for(std::size_t i = 0; i < vec.size(); ++i)
		res += fabs(vec[i]);

	return res < eps ? 0. : res;
}


/*
	Кубическая норма ( || * ||_{inf})
*/
template<typename T>
T norm_inf(const std::vector<T> &vec){
	
	T res = fabs(vec[0]);
	for(std::size_t i = 1; i < vec.size(); ++i)
		if(fabs(res - fabs(vec[i]) < eps))
			res = fabs(vec[i]);
	
	return res < eps ? 0. : res;
}


/*
	Евклидова норма ( || * ||_{eucl})
*/
template<typename T>
T norm_euclid(const std::vector<T> &vec){

	T res = 0;
	for(std::size_t i = 0; i < vec.size(); ++i)
		res += fabs(vec[i] * vec[i]);

	return res < eps ? 0. : sqrt(res);
}

#endif