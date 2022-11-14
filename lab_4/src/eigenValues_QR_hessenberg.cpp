
#include "eigenValues_QR_hessenberg.h"

#include "Matrix_public.h"

#include <iostream>
#include "QR_method.h"



double* hessenberg(double**& A, size_t& n){

    double g = 0.0;

    double* lambda_mas = new double [n];

    double** Ak = new double* [n];
    for (size_t j = 0; j < n; ++j){
        Ak[j] = new double[n];
    }
    copy_Matrix_A_to_B(A, Ak, n);

    double a, b, c, s;

    print_Matrix(Ak, n);

	for (size_t i = 1; i < n-1; ++i) {

		for (size_t j = i + 1; j < n; ++j) {

			coefs_hessen(i, j, c, s, Ak, n);

			for (size_t k = 0; k < n; ++k) {
                
				a = Ak[i][k];
				b = Ak[j][k];

				Ak[i][k] = (c * a + s * b);
				Ak[j][k] = (-s * a + c * b);

                // a = A[i][k];
				// b = A[j][k];

				// if (a == 0 && b != 0){
				// 	Ak[i][k] = (s * b);
				// 	Ak[j][k] = (c * b);
				// }
				// if (a != 0 && b == 0){
				// 	Ak[i][k] = (c * a);
				// 	Ak[j][k] = (-s * a );
				// }

				// if (a == 0 && b == 0){
				// 	Ak[i][k] = 0;
				// 	Ak[j][k] = 0;
				// }

				// if (a != 0 && b != 0){
					
				// 	Ak[i][k] = (c * a + s * b);
				// 	Ak[j][k] = (-s * a + c * b);

				// }
				
			}	
            for (size_t k = 0; k < n; ++k) {
                

                
                a = Ak[k][i];
				b = Ak[k][j];

				Ak[k][i] = (c * a + s * b);
				Ak[k][j] = (-s * a + c * b);

                // a = A[i][k];
				// // b = A[j][k];

				// if (a == 0 && b != 0){
				// 	Ak[k][i] = (s * b);
				// 	Ak[k][j] = (c * b);
				// }
				// if (a != 0 && b == 0){
				// 	Ak[k][i] = (c * a);
				// 	Ak[k][j] = (-s * a );
				// }

				// if (a == 0 && b == 0){
				// 	Ak[k][i] = 0;
				// 	Ak[k][j] = 0;
				// }

				// if (a != 0 && b != 0){
					
				// 	Ak[k][i] = (c * a + s * b);
				// 	Ak[k][j] = (-s * a + c * b);

				// }



            }

    		//print_Matrix(Ak, n);

        
		}


	}

    check_matrix_zero(Ak, n);
    print_Matrix(Ak, n);


    double** Q = new double* [n];
    double** R = new double* [n];
    //double** Eg = new double* [n];



    for (size_t j = 0; j < n; ++j){
        //Ak[j] = new double[n];

        Q[j] = new double[n];
        R[j] = new double[n];
    }

    

    //copy_Matrix_A_to_B(A, Ak, n);

    int iter = 0;
    double sum = 0;

    while (true)
	//for(size_t u = 0; u < 5; ++u)
    {   
        ++iter;
        g = Ak[n-1][n-1];

        for (size_t i = 0; i < n; ++i){
            Ak[i][i] -= g;
        }

        QR_matrix(Q,R,Ak,n);

        matrix_mult(R, Q, Ak, n);
        
        for (size_t i = 0; i < n; ++i){
            Ak[i][i] += g;
        }

        sum = 0;

        for (size_t i = 0; i < n-2; ++i) {
			sum += fabs(Ak[i+1][i]);               
        }

        if (sum < eps)
            break;

		//print_Matrix(Ak, n);


    } 
    
    print_Matrix(Ak, n);

    for(size_t i = 0; i < n; ++i){
        lambda_mas[i] = Ak[i][i];
        //cout << "Lambda_" << i << " = " << Ak[i][i] << "\n";
    }


    cout << "iter = " << iter << "\n";

    return lambda_mas;
}



int coefs_hessen(size_t& k, size_t& l, double& alpha, double& beta, double**& A, size_t& n) {

	if (k < n && l < n) {

		double temp = sqrt(pow(A[k][k-1], 2) + pow(A[l][k-1], 2));

		alpha = A[k][k-1] / temp;

		beta = A[l][k-1] / temp;

		if (fabs(alpha) < eps)
			alpha = (1 / alpha > 0) ? alpha : alpha*(-1);
		if (fabs(beta) < eps)
			beta = (1 / beta > 0) ? beta : beta*(-1);

	}


	return 0;
}