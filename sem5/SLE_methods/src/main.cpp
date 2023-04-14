#include <iostream>

#include "../../../structures/linalg/Matrix_n_Vector.hpp"
#include "./include/Solver_SLE.hpp"

int main(int args, char** argv) {
    std::string path = "../test_files/P_DATA6.TXT";

    Solver_SLE<double> solver;

    // Matrix<double> K;
    // Vector<double> r;

    // read_System(path, K, r);

    // std::cout << K << "\n";
    // std::cout << r << "\n\n";

    // std::cout << solver.Gauss(K, r);

    // auto xQR = solver.QR(K, r);

    // std::cout << "QR\n\n";
    // std::cout << std::get<0>(xQR) << "\n";
    // std::cout << std::get<1>(xQR) << "\n";
    // std::cout << std::get<2>(xQR) << "\n";


    //Banded matrix test

    /*   Ax = b

        Regular matrix A:
        [b, c, 0,    ....]   [d]
        [a, b, c,    ....]   [d]
        [0, a, b, c, ....]   [d]
        .................. = ...
        [......., a, b, c]   [d]
        [.........., a, b]   [d]  

        b = [d, d, d, d, d]
        Banded matrix from A:

        banded_matrix = [
            [c, c, c, c, *],
            [b, b, b, b, b],
            [*, a, a, a, a],
        ]
        
        */

    Matrix<double> A(3, 4);
    Vector<double> b(4);
    A[0][0] = 2;
    A[0][1] = 4;
    A[0][2] = 5;

    A[1][0] = 1;
    A[1][1] = 3;
    A[1][2] = 4;
    A[1][3] = 9;

    
    A[2][1] = 2;
    A[2][2] = 6;
    A[2][3] = 8;

    A.print();

    b[0] = 1;
    b[1] = 3;
    b[2] = 4;
    b[3] = 5;

    auto x = solver.solve_banded(1, 1, A, b);

    std::cout << x << "\n";

    return 0;
}
