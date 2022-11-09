#include "inverse_iter_method.h"
#include "Matrix_public.h"
#include "QR_method.h"
#include "method_gauss.h"
#include "Relaxation_method.h"

#include <iostream>
#include <cmath>
#define eps 1e-15


using namespace std;



double* inverse_iter_method(double**& A, double*& x0, double lambda, size_t n){


    double* x = new double[n];
    double* xk = new double[n];
    double* tmp = new double [n];

    double** C = new double*[n];

    for (size_t j = 0; j < n; ++j){
			C[j] = new double[n + 1];
    }

    copy_Sys_A_to_B(A, C, n);

    for (size_t i = 0; i < n; ++i){
        C[i][i] -= lambda;
    }

    double norm = 0.0;   
    
    copy_Vector_A_to_B(x0, xk, n);
    copy_Vector_A_to_B(x0, x, n);


    double** Q = new double* [n];
    double** R = new double* [n];

    for (size_t j = 0; j < n; ++j){
        Q[j] = new double[n];
        R[j] = new double[n];
    }

    QR_matrix(Q, R, C, n);

    copy_Matrix_A_to_B(R, C, n);
    transpose_matrix(Q, n);

    check_sys_matrix_zero(C, n);
    
    int iter = 0;
    double sum = 0;
    
    while (true)
    {
        ++iter;
        
        matrix_vector_mult(Q, x, tmp, n);

        for(size_t i = 0; i < n; ++i){
            C[i][n] = tmp[i];
        }

        check_sys_matrix_zero(C, n);
        reverse_course(C, n, xk);
        
        e_norm(xk, norm, n);

        for(size_t i = 0; i < n; ++i){
            xk[i] /= norm;
        }

        sum = 0;
        for(size_t i = 0; i < n; ++i){
            tmp[i] =  ((x[i]) - (xk[i]));
            
            
            sum += fabs(x[i] + xk[i]);
            
        }  
        copy_Vector_A_to_B(xk, x, n);

        if(sum < eps){
            break;
        }

        norm_inf_vector(tmp, norm, n);
        if (norm < eps){
            break;
        }

    }



    cout << "k = " << iter << "\n";

    for (size_t j = 0; j < n; ++j){
        delete[] C[j];
        delete[] Q[j];
        delete[] R[j];


    }
    delete[] C;
    delete[] Q;
    delete[] R;
    delete[] xk;
    delete[] tmp;

    

    return x;
}


double* inverse_iter_method_Rayleigh(double**& A, double*& x0, double lambda, size_t n){

    double* x = new double[n];
    double* xk = new double[n];
    double* tmp = new double [n];

    double** C = new double*[n];

    for (size_t j = 0; j < n; ++j){
			C[j] = new double[n + 1];
    }
    copy_Sys_A_to_B(A, C, n);

    double norm = 0.0;   
    
    copy_Vector_A_to_B(x0, xk, n);
    copy_Vector_A_to_B(x0, x, n);


    double sum = 0;
    int iter = 0;

    while (true)
    {
        ++iter;

        for (size_t i = 0; i < n; ++i){
            C[i][i] -= lambda;
            C[i][n] = x[i];

        }

        if (false == inv_check(C, n)){
            return x;
        }
        
        method_gauss(C, n, xk);

        e_norm(xk, norm, n);


        for(size_t i = 0; i < n; ++i){
            xk[i] /= norm;
        }
        sum = 0;
        for(size_t i = 0; i < n; ++i){
            tmp[i] = x[i] - xk[i];
            sum += fabs(x[i] + xk[i]);

        }    

        copy_Vector_A_to_B(xk, x, n);

        if(sum < eps){
            break;
        }

        norm_inf_vector(tmp, norm, n);

        if (norm < eps){
            break;
        }

        matrix_vector_mult(C, xk, tmp, n);
        scalar_mult(tmp, xk, lambda, n);

    }

    cout << "k = " << iter << "\n";

    for (size_t j = 0; j < n; ++j){
        delete[] C[j];
    }
    delete[] C;
    delete[] xk;
    delete[] tmp;

    return x;
}