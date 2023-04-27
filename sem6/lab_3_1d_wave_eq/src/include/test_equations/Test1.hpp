#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_1d_wave_equation.hpp"

/**
 * Тест 1
*/
template<class T>
class Test1: virtual public Class_1d_wave_equation<T>{

    public:

        Test1(){
            
            //test name
            this->_name = std::string (__func__);

            this->DEFAULT_TEST();
        }

};