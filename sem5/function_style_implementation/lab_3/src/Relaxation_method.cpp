
#include <iostream>
#include "Matrix_public.h"
#include <cmath>

double *Relaxation_method(double **&A, double *&x0, double &w, size_t n)
{

    // double* xt = new double[n];
    // xt[0] = 5;
    // xt[1] = -7;
    // xt[2] = 12;
    // xt[3] = 4;


    double *x = new double[n];
    double *xk = new double[n];
    double *tmp = new double[n];

    double **C = new double *[n];

    double **L = new double *[n];
    double **D = new double *[n];
    double **U = new double *[n];

    for (size_t i = 0; i < n; ++i)
    {
        L[i] = new double[n];
        D[i] = new double[n];
        U[i] = new double[n];
    }
    LDU_method(A, L, D, U, n);

    double **TMP1 = new double *[n];
    double **TMP2 = new double *[n];

    for (size_t j = 0; j < n; ++j)
    {
        C[j] = new double[n + 1];
        TMP1[j] = new double[n];
        TMP2[j] = new double[n];
    }

    copy_Sys_A_to_B(A, C, n);

    matrix_make_E(TMP1, n);

    copy_Sys_A_to_B(A, TMP1, n);

    // cout << " -------------- \n";
    // print_Matrix(TMP1, n);
    // //print_Matrix(TMP2, n);
    // cout << " -------------- \n";

    matrix_number_mult(L, w, n);
    matrix_number_mult(TMP1, w, n);

    // cout << " -------------- \n";
    // print_Matrix(L, n);
    // print_Matrix(TMP1, n);
    // cout << " -------------- \n";

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {

            C[i][j] = D[i][j] + L[i][j];
            TMP2[i][j] = D[i][j] + L[i][j] - TMP1[i][j];
        }
    }



    reverse_matrix_qr(C, TMP1, n);


    matrix_mult(TMP1, TMP2, C, n);



    //print_Matrix(C, n);
    LDU_method(C, L, D, U, n);

    double dC = 0.0;
    norm_inf_matrix(C, dC, n);
    double normC = dC;

    //cout << "|C| = " << dC << "\n";


    double dCU = 0.0;
    norm_inf_matrix(U, dCU, n);

    dC = (1 - dC) / dCU;

    double sum = 0.0;
    double norm = 0.0;

    copy_Vector_A_to_B(x0, xk, n);

    int iter = 0;
    double normX = 0.0;   

    //for(int p = 0; p < 61; ++p)
    while (true)
    {
        ++iter;
        for (size_t i = 0; i < n; ++i)
        {

            
                sum = 0.0;
                if (i > 0)
                    for (size_t j = 0; j <= i - 1; ++j)
                    {
                        sum -= A[i][j] * x[j];
                    }
                if (i < n)
                    for (size_t j = i+1; j < n; ++j)
                    {
                        sum -= A[i][j] * xk[j];
                    }
                x[i] = (1 - w) * xk[i] + w * (sum + A[i][n]) / A[i][i];

        }

        for (size_t i = 0; i < n; ++i)
        {
            tmp[i] = x[i] - xk[i];
        }
        norm_inf_vector(tmp, norm, n);



        // matrix_vector_mult(A, x, tmp, n);
        // for(size_t i = 0; i < n; ++i)
        //     tmp[i] -= A[i][n];



        // if (norm < eps)
        // {
        //     break;
        // }

        // if (norm < dC * eps)
        // {
        //     break;
        // }

    
        norm_inf_vector(x, normX, n);

        if (norm / (normX + 0.01) < eps)
        {
            break;
        }
        

        copy_Vector_A_to_B(x, xk, n);

        // if (iter == 1){
        //     cout << "k_est = " << ceil(log((1 - normC) * eps / norm) / log(normC) )<< endl;
        // }

    }

    // cout << "k = " << iter << endl;

    // for(size_t i = 0; i < n; ++i){
    //         tmp[i] = x[i] - xt[i];
    //     }              
    // norm_inf_vector(tmp, norm, n);

    // cout << "|norm_err| = " << norm << "\n";

    for (size_t j = 0; j < n; ++j)
    {
        delete[] C[j];
        // delete[] TMP1[j];
        // delete[] TMP2[j];
        delete[] L[j];
        delete[] D[j];
        delete[] U[j];

    }

    delete[] C;
    delete[] TMP1;
    delete[] TMP2;
    delete[] L;
    delete[] D;
    delete[] U;


    delete[] xk;
    delete[] tmp;

    return x;
}



#include <iostream>
#include "Matrix_public.h"

double *Relaxation_method(double*& a, double*& b, double*& c, double*& d, double*& x0, double& w, size_t n)
{

    double *x = new double[n];
    double *xk = new double[n];
    double *tmp = new double[n];

    double **C = new double *[n];

    double **L = new double *[n];
    double **D = new double *[n];
    double **U = new double *[n];

    for (size_t i = 0; i < n; ++i)
    {
        L[i] = new double[n];
        D[i] = new double[n];
        U[i] = new double[n];
    }
    //LDU_method(A, L, D, U, n);

    double **TMP1 = new double *[n];
    double **TMP2 = new double *[n];

    for (size_t j = 0; j < n; ++j)
    {
        C[j] = new double[n + 1];
        TMP1[j] = new double[n];
        TMP2[j] = new double[n];
    }

    matrix_make_E(C, n);

    for(size_t i = 0; i < n; ++i){
        
        if(i > 0){
            C[i][i-1] = a[i];
        }
        C[i][i] = b[i];
        if(i < n - 1){
            C[i][i+1] = c[i];
        }        
        C[i][n] = d[i];

    }
    //print_Sys_Matrix(C, n);
    LDU_method(C, L, D, U, n);

    //copy_Sys_A_to_B(A, C, n);

    matrix_make_E(TMP1, n);

    copy_Sys_A_to_B(C, TMP1, n);

    // cout << " -------------- \n";
    // print_Matrix(TMP1, n);
    // //print_Matrix(TMP2, n);
    // cout << " -------------- \n";

    matrix_number_mult(L, w, n);
    matrix_number_mult(TMP1, w, n);

    // cout << " -------------- \n";
    // print_Matrix(L, n);
    // print_Matrix(TMP1, n);
    // cout << " -------------- \n";

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {

            C[i][j] = D[i][j] + L[i][j];
            TMP2[i][j] = D[i][j] + L[i][j] - TMP1[i][j];
        }
    }

    reverse_matrix_qr(C, TMP1, n);
    //print_Matrix(C, n);

    matrix_mult(TMP1, TMP2, C, n);

    //print_Matrix(C, n);

    LDU_method(C, L, D, U, n);
    //print_Matrix(U, n);
    double dC = 0.0;
    norm_inf_matrix(C, dC, n);
    //cout << "|C| = " << dC << "\n";
    double normC = dC;


    double dCU = 0.0;
    norm_inf_matrix(U, dCU, n);

    dC = (1 - dC) / dCU;
    //cout << "Парам останова = " << dC << "\n";


    double sum = 0.0;
    double norm = 0.0;

    copy_Vector_A_to_B(x0, xk, n);
    int iter = 0;
    while (true)
    {
        ++iter;
        for (size_t i = 0; i < n; ++i)
        {
            sum = 0.0;
            if (i > 0)
            {
                sum -= a[i] * x[i-1];
            }
            if (i < n-1)
            {
                sum -= c[i] * xk[i+1];
            }
            x[i] = (1 - w) * xk[i] + w * (sum + d[i]) / b[i];

        }

        for (size_t i = 0; i < n; ++i)
        {
            tmp[i] = x[i] - xk[i];
        }

        norm_inf_vector(tmp, norm, n);

        if (norm < dC * eps)
        {
            break;
        }
        copy_Vector_A_to_B(x, xk, n);

        // if (iter == 1){
        //     cout << "k_est = " << ceil(log((1 - normC) * eps / norm) / log(normC) )<< endl;
        // }

    }

    //cout << "k = " << iter << endl;

    for (size_t j = 0; j < n; ++j)
    {
        delete[] C[j];
        delete[] TMP1[j];
        delete[] TMP2[j];
        delete[] L[j];
        delete[] D[j];
        delete[] U[j];

    }

    delete[] C;
    delete[] TMP1;
    delete[] TMP2;
    delete[] L;
    delete[] D;
    delete[] U;


    delete[] xk;
    delete[] tmp;

    return x;
}