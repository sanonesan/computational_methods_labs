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

        Test1<T> set_01(){
            
            //test name
            this->_name = "Test1_01";
            this->DEFAULT_TEST((T)0.1);

            return *this;
        }

};