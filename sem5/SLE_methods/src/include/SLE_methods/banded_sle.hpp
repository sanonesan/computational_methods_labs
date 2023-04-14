#ifndef BANDED_SLE_HPP
#define BANDED_SLE_HPP

#include <iostream>
#include <cmath>

#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"

template<typename T>
Vector<T> banded_sle(const std::size_t l, const std::size_t u, Matrix<T>& banded_matrix, Vector<T>& b){

    Vector<T> solution(b);
    
    if(u == 1 && l == 1){
        /*

        banded_matrix = [
            [c, c, c, c, *],
            [b, b, b, b, b],
            [*, a, a, a, a],
        ]

        [b, c, 0,    ....]   [d]
        [a, b, c,    ....]   [d]
        [0, a, b, c, ....]   [d]
        .................. = ...
        [......., a, b, c]   [d]
        [.........., a, b]   [d]        
        
        
        */
        std::vector<T> alpha;
        std::vector<T> beta;
        T tmp = 0;

        alpha.push_back( - banded_matrix[0][0] / banded_matrix[1][0] );
        beta.push_back(b[0] / banded_matrix[1][0]);

        for(std::size_t i = 1; i < b.size(); ++i){
            tmp = banded_matrix[1][i] + banded_matrix[2][i] * alpha[i - 1];
            (i < b.size() - 1) ? alpha.push_back( - banded_matrix[0][i] / tmp ) : alpha.push_back(0.);
            beta.push_back( (b[i] - banded_matrix[2][i] * beta[i-1]) / tmp );
        }
        
        solution[b.size() - 1] = beta[b.size() - 1];
        
        for(std::size_t j = b.size() - 2; j >= 0; --j){
            solution[j] = beta[j] + alpha[j] * solution[j + 1];
            if(j == 0){
                break;
            }
        }
        
    }  
    
    return solution;

};


#endif