
#include <iostream>
#include "Matrix_public.h"


double* simple_iter_method(double**& A, double*& x0, double tau, size_t n){


    double* x = new double[n];
    double* xk = new double[n];
    double* tmp = new double [n];

    double** C = new double*[n];
    for (size_t j = 0; j < n; ++j)
			C[j] = new double[n + 1];
    
    copy_Sys_A_to_B(A, C, n);

    for (size_t i = 0; i < n; ++i){
        for (size_t j = 0; j < n + 1; ++j){
            C[i][j] *= (j == n) ? tau : -tau;           
        }
        C[i][i] += 1;
    }

    //print_Sys_Matrix(C, n);

    double dC = 0.0;
    norm_inf_matrix(C, dC, n);
    cout << "|C| = " << dC << "\n";

    dC = (1 - dC) / dC;
    
    copy_Vector_A_to_B(x0, xk, n);
    copy_Vector_A_to_B(x0, x, n);

    double norm = 0.0;   
    int iter = 0;

    while (true)    
    {
        ++iter;
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

    cout << "k = " << iter << "\n";

    for (size_t j = 0; j < n; ++j){
        delete [] C[j];
    }
    delete [] C;
    delete [] xk;
    delete [] tmp;
    

    return x;
}