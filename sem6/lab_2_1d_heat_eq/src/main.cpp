
#include<iostream>
#include"../../../structures/linalg/Matrix_n_Vector.hpp"
int main(int args, char **argv){

	
    Vector<double> v({1., 4., 5.});
    Vector<double> w({1., 3., 5.});
    Vector<double> k(w);
    k += v;
    //w += 3;
    // // std::cout << v * w  << "\t" << k  << "\t" << w;

    typedef double T;

	Matrix<T> K, M;

	K.read_Square_Matrix("../../../sem5/lab_1/test_files/test.txt");
	M.read_Square_Matrix("../../../sem5/lab_1/test_files/test.txt");

	// K._array[0][0] = 1;	
	// K._array[0][1] = 2;	
	// K._array[1][0] = 1;	
	// K._array[1][1] = 2;
	// K._array[2][0] = 1;	
	// K._array[2][1] = 2;

	// Vector<T> m(2);
	std::cout << K << std::endl;
	std::cout << M << std::endl;

    K * M;
	std::cout << K << std::endl;

    K *= M;
    K - M;
    K += M;
    K -= M;

    
	std::cout << (K.dot(w)) << std::endl;
    

	// m[0] = 1;
	// m[1] = 2;

	// std::cout << m.dot(K) << std::endl;

	Matrix<T> Z(3,3);
	std::cout << Z[0][0] << std::endl;

	
    return 0;
}
