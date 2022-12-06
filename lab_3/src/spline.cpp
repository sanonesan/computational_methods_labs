#include "spline.h"
#include "Matrix_public.h"
#include "Relaxation_method.h"
#include "method_gauss.h"

using namespace std;

double** spline(double *&x_mas, double *&y_mas, size_t &n)
{

    double spline = 0;

    double a, b, c, d;
    size_t k = n - 1;
    double **A = new double *[k];
    for (size_t i = 0; i < k; ++i)
    {
        A[i] = new double[4];
    }

    double *h = new double[k];
    double *g = new double[k];

    for (size_t i = 0; i < k; ++i)
    {
        h[i] = x_mas[i+1] - x_mas[i];

        g[i] = (y_mas[i + 1] - y_mas[i]) / h[i];
    }

    size_t kk = k -1;

    double *ai = new double[kk];
    double *bi = new double[kk];
    double *ci = new double[kk];
    double *di = new double[kk];

    ai[0] = 0;
    ci[kk-1] = 0;


    double** C = new double*[kk];
    for(size_t i = 0; i < k; ++i){
        C[i] = new double [kk + 1];
    }

    // for(size_t i = 0; i < kk; ++i){
    //     for(size_t j = 0; j < kk; ++j){
            
    //         if(i == j){
    //             C[i][j] = 2 * (h[i] + h[i+1]); //bi
    //             bi[i] = 2 * (h[i] + h[i+1]);
    //         }
    //         else {
    //             if(i - 1 == j && i - 1 >= 0){
    //                 C[i][j] = h[i - 1]; //ai
    //                 ai[i] = h[i-1];
    //             }
    //             else {
    //                 if(i + 1 == j && i + 1 <= n){
    //                     C[i][j] = h[i + 1]; //ci
    //                     ci[i] = h[i+1];
    //                 }
    //                 else{
    //                     C[i][j] = 0.0;
    //                 }
    //             }
    //         }
    //     }
    //         C[i][kk] = 3 * (g[i+1] - g[i]);
    //         di[i] = 3 * (g[i+1] - g[i]);
    // }

    for(size_t i = 0; i < kk; ++i){
        bi[i] = 2 * (h[i] + h[i+1]);
        ai[i] = h[i-1];
        ci[i] = h[i+1];
        di[i] = 3 * (g[i+1] - g[i]);
    }


    double* c_res1 = new double[kk];
    double w = 1.0;
    double* x0 = new double [kk];
	for (size_t i = 0; i < kk; ++i){
		x0[i] = 0;
	}

    //c_res1 = Relaxation_method(C, x0, w, kk);
    c_res1 = Relaxation_method(ai, bi, ci, di, x0, w, kk);

    //print_Vector(c_res1, kk);
    //method_gauss(C, kk, c_res1);

    double* c_res = new double[n];

    for(size_t i = 0; i < k; ++i){
        c_res[i] = (i == 0) ? 0.0 : c_res1[i-1];
    }
    // print_Vector(c_res1, kk);
    // print_Vector(c_res, n);

    delete[] c_res1;

    for(size_t i = 0; i < k; ++i){
        A[i][0] = y_mas[i];
        A[i][1] = g[i] - (c_res[i+1] + 2 * c_res[i]) * h[i] / 3;
        A[i][2] = c_res[i];
        A[i][3] = (c_res[i+1] - c_res[i]) / 3 / h[i];
    }

    // for(size_t i = 0; i < k; ++i){
    //     for(size_t j = 0; j < 4; ++j){
    //         cout << A[i][j] << "    ";
    //     }
    //     cout << "\n";
    // }

    return A;
}