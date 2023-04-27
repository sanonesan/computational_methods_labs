#include <iostream>

#include "../../../structures/linalg/Matrix_n_Vector.hpp"
#include "./include/Solver_SLE.hpp"

int main(int args, char** argv) {

    std::string path = "../test_files/input12.txt";

    Solver_SLE<double> solver;

    Matrix<double> K;
    Vector<double> r;

    read_System(path, K, r);

    // std::cout << K << "\n";
    // std::cout << r << "\n\n";

    std::cout << solver.Gauss(K, r) << "\n";

    auto xQR = solver.QR(K, r);

    // std::cout << "QR\n\n";
    std::cout << std::get<0>(xQR) << "\n";
    // std::cout << std::get<1>(xQR) << "\n";
    // std::cout << std::get<2>(xQR) << "\n";
    
    Vector<double> x0(4);
    std::cout << solver.solve_Jacoby(K, r, x0) << "\n";
    std::cout << solver.solve_simple_iter(K, r, x0, 0.005, 10);
    
    

    return 0;
}
