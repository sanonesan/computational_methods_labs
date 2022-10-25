
#include <iostream>
#include "Matrix_public.h"


double* Jacoby_method(double**& A, double*& x0, size_t n){


    double* x = new double[n];
    double* xk = new double[n];
    double* tmp = new double [n];

    double** C = new double*[n];
    for (size_t j = 0; j < n; ++j){
			C[j] = new double[n + 1];
    }
    copy_Sys_A_to_B(A, C, n);    

    for (size_t i = 0; i < n; ++i){
        for (size_t j = 0; j < n+1; ++ j){
            C[i][j] *= (i == j && j < n) ? 0.0 : - 1 / A[i][i];
        }
        C[i][n] = A[i][n] / A[i][i];
    }
    //print_Sys_Matrix(C, n);

    for (size_t i = 0; i < n; ++i){
        
        x[i] = A[i][n];
    }   

    double dC = 0.0;
    norm_inf_matrix(C, dC, n);
    cout << "|C| = " << dC << "\n";

    dC = (1 - dC) / dC;

    copy_Vector_A_to_B(x0, xk, n);
    copy_Vector_A_to_B(x0, x, n);

    double norm = 0.0;   

    while (true)    
    {
        matrix_vector_mult(C, xk, tmp, n);
        for(size_t i = 0; i < n; ++i){
            x[i] = tmp[i] + C[i][n];
        } 
        
        for(size_t i = 0; i < n; ++i){
            tmp[i] = x[i] - xk[i];
        }              
        norm_inf_vector(tmp, norm, n);

        if (norm < dC * eps){
            break;
        }

        copy_Vector_A_to_B(x, xk, n);

    }

    for (size_t j = 0; j < n; ++j){
        delete [] C[j];
    }
    delete [] C;
    delete [] xk;
    delete [] tmp;   

    return x;
}