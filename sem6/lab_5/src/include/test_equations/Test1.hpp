#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_Fred.hpp"

/**
 * Тест 1
*/
template<class T>
class Test1: virtual public Class_Fred<T>{

    public:

        Test1(){
            
            //test name
            this->_name = std::string (__func__);

            this->DEFAULT_TEST();
        }

};