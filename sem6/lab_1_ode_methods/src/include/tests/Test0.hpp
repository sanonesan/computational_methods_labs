#pragma once

#include<vector>
#include<cmath>
#include<functional> 

/*
// ...........Тест..из..методички.............. //
*/
template<typename T>
class Test0{

    public:

        std::vector<T> _x0 = {0.5};
        std::vector<std::function<T (const std::vector<T>& x, const T t)>> _ode_system;


        Test0(){
            
            // dx/dt
            auto func_1 = [](const std::vector<T>& x, const T t) -> T{
                return cos(t) - x[0];
            };


            _ode_system.push_back(func_1);
            // _ode_system.push_back(func_2);
        }

};
