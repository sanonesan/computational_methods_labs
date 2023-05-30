#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_2d_Poisson_equation_in_rectangle.hpp"

/**
 * Тест 1
*/
template<class T>
class Test1: virtual public Class_2d_Poisson_equation_in_rectangle<T>{

    public:

        Test1(){
            
            //test name
            this->_name = std::string (__func__);

            this->DEFAULT_TEST();
        }

};