#include "interpol.h"

double interpol_lagrange(double x, double *&x_mas, double *&y_mas, size_t &n)
{
    double lagrange_pol = 0;
    double c_i = 0;

    for (int i = 0; i < n; i++)
    {
        c_i = 1;
        for (int j = 0; j < n; j++)
        {
            if (j == i)
                continue;
            c_i *= (x - x_mas[j]) / (x_mas[i] - x_mas[j]);
        }
        lagrange_pol += c_i * y_mas[i];
    }
    return lagrange_pol;
}
