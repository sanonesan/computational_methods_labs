#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_Fred.hpp"

/**
 * Тест 1
*/
template<class T>
class Test2: virtual public Class_Fred<T>{

    public:

        Test2(){
            
            //test name
            this->_name = std::string (__func__);

            this->DEFAULT_TEST_1();
        }

        Test2<T> set_01(){
            
            //test name
            this->_name = "Test2_01";

            this->DEFAULT_TEST_1((T)0.1);

            return *this;
        }

};