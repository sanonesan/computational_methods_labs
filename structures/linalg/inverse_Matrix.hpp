// #include "Matrix_n_Vector.hpp"

// #include <tuple>

// #include "../../sem5/lab_1_SLE_methods/src/include/Solver_SLE.hpp"
// #include "../../sem5/lab_1_SLE_methods/src/include/SLE_methods/method_gauss.hpp"

// template <class T>
// Matrix<T> inv(Matrix<T> A) {


//     // Matrix<T> A(*this);
//     Solver_SLE<T> solver;
//     Matrix<T> InvMatrix(A);

//     auto QR_answer = solver.QR_decomposion(A);

//     std::size_t n = A.get_rows();
//     Vector<T> b(n);
//     Vector<T> tmp(n);
//     Vector<T> tmp_solution(n);
//     Matrix<T> B(std::get<0>(QR_answer));

    
//     for (std::size_t i = 0; i < n; ++i){
        
//         for (std::size_t j = 0; j < n; ++j){
//             b[j] = (i == j) ? 1. : 0. ;
//         }

//         tmp = B.dot(b);

//         reverse_course(std::get<1>(QR_answer), tmp, InvMatrix[i]);

//     }
//     InvMatrix.transpose_this();
//     InvMatrix.check_matrix_zero();
    
//     return InvMatrix;
// }