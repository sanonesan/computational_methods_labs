
#include <iostream>
#include "Matrix_public.h"
#include <cmath>

double* simple_iter_method(double**& A, double*& x0, double tau, size_t n){

    double* xt = new double[n];
    xt[0] = 5;
    xt[1] = -7;
    xt[2] = 12;
    xt[3] = 4;


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
        tmp[i] = C[i][n];
    }

    //print_Sys_Matrix(C, n);
    
    double dC = 0.0;
    double normC = 0.0;
    double norm = 0.0;   

    norm_inf_matrix(C, dC, n);
    normC = dC;
    cout << "|C| = " << dC << "\n";


    //norm_inf_vector(tmp, norm, n);
    //cout << "|tau*b| = " << norm << "\n";

    dC = (1 - dC) / dC;
    
    copy_Vector_A_to_B(x0, xk, n);
    copy_Vector_A_to_B(x0, x, n);

    int iter = 0;
    double normX = 0.0;   

    //for(int p = 0; p < 604; ++p)
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

        // if ((norm - dC * eps) < eps){
        //     break;
        // }

        // matrix_vector_mult(A, x, tmp, n);
        // for(size_t i = 0; i < n; ++i)
        //     tmp[i] -= A[i][n];

        // norm_inf_vector(tmp, norm, n);

        // if (norm < eps)
        // {
        //     break;
        // }

        norm_inf_vector(x, normX, n);

        if (norm / (normX + 0.01) < eps)
        {
            break;
        }

        copy_Vector_A_to_B(x, xk, n);
        
        if (iter == 1){
            cout << "k_est = " << ceil(log((1 - normC) * eps / norm) / log(normC) )<< endl;
        }
    }

    cout << "k = " << iter << "\n";
    

    for(size_t i = 0; i < n; ++i){
            tmp[i] = x[i] - xt[i];
        }              
    norm_inf_vector(tmp, norm, n);

    cout << "|norm_err| = " << norm << "\n";

    for (size_t j = 0; j < n; ++j){
        delete [] C[j];
    }
    delete [] C;
    delete [] xk;
    delete [] tmp;
    

    return x;
}