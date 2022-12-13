#include "spline.h"
#include "Matrix_public.h"
#include "Relaxation_method.h"
#include "method_gauss.h"

using namespace std;

double** spline(double *&x_mas, double *&y_mas, size_t &n)
{
    double spline = 0;

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

    for(size_t i = 0; i < kk; ++i){
        bi[i] = 2 * (h[i] + h[i+1]);
        ai[i] = h[i-1];
        ci[i] = h[i+1];
        di[i] = 3 * (g[i+1] - g[i]);
    }

    double* c_res1 = new double[kk];
    
    triagonal_matrix_algorithm(ai, bi, ci, di, c_res1, kk);

    double* c_res = new double[n];

    for(size_t i = 0; i < k; ++i){
        c_res[i] = (i == 0) ? 0.0 : c_res1[i-1];
    }

    for(size_t i = 0; i < k; ++i){
        A[i][0] = y_mas[i];
        A[i][1] = g[i] - (c_res[i+1] + 2 * c_res[i]) * h[i] / 3;
        A[i][2] = c_res[i];
        A[i][3] = (c_res[i+1] - c_res[i]) / 3 / h[i];
    }

    delete[] ai;
    delete[] bi;
    delete[] ci;
    delete[] di;

    delete[] g;
    delete[] h;
    delete[] c_res;
    delete[] c_res1;
    
    return A;
}