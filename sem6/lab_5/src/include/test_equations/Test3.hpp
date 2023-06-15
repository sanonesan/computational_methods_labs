#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_Fred.hpp"

/**
 * Тест 1
*/
template<class T>
class Test3: virtual public Class_Fred<T>{

    public:

        Test3(int N = 10){
            
            //test name
            this->_name = std::string (__func__);
            this->_name += "_" + std::to_string(N);
            this->DEFAULT_TEST_singular(N);
        }

};