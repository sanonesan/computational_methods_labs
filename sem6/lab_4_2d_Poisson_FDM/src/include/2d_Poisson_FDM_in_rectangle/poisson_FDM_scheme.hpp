#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>

#include "../../../../../sem5/SLE_methods/src/include/Solver_SLE.hpp"
#include "../../../../../structures/linalg/Matrix_n_Vector.hpp"
#include "../Class_2d_Poisson_equation_in_rectangle.hpp"

template <typename T>
void poisson_FDM_scheme(

    const Class_2d_Poisson_equation_in_rectangle<T>& wave_equation,
    const std::vector<T>& time,
    const std::vector<T>& x,
    std::vector<T>& y_j_minus_1,
    const T& tolerance,
    std::ofstream &fout){
    
    
};
