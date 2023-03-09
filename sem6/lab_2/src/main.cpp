
#include<iostream>
#include"../../../structures/linalg/Matrix_n_Vector.hpp"


int main(int args, char **argv){

	typedef double T;

	Matrix<T> K(2), M;
	//K.read_Square_Matrix("../../../sem5/lab_1/test_files/test.txt");
	//M.read_Square_Matrix("../../../sem5/lab_1/test_files/test.txt");

	K._array[0][0] = 1;	
	K._array[0][1] = 2;	
	K._array[1][0] = 1;	
	K._array[1][1] = 2;
	// K._array[2][0] = 1;	
	// K._array[2][1] = 2;

	Vector<T> m(2);

	m._vector[0] = 1;
	m._vector[1] = 2;

	std::cout << K.dot(m) << std::endl;
	std::cout << m.dot(K) << std::endl;

	
    return 0;
}

