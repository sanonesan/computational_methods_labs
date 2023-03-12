#include <iostream>

#include "../../../structures/linalg/Matrix_n_Vector.hpp"
#include "./include/Solver_SLE.hpp"

int main(int args, char** argv) {
    std::string path = "../test_files/P_DATA6.TXT";

    Solver_SLE<double> solver;

    Matrix<double> K;
    Vector<double> r;

    read_System(path, K, r);

    std::cout << K << "\n";
    std::cout << r << "\n\n";

    std::cout << solver.Gauss(K, r);

    auto xQR = solver.QR(K, r);

    std::cout << "QR\n\n";
    std::cout << std::get<0>(xQR) << "\n";
    std::cout << std::get<1>(xQR) << "\n";
    std::cout << std::get<2>(xQR) << "\n";

    return 0;
}
