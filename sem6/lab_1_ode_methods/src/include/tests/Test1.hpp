#pragma once

#include<vector>
#include<cmath>
#include<functional> 

/*
// ...........Тест..из..методички.............. //
*/
template<typename T>
class Test1{

    public:

        std::vector<T> _x0 = {1.,0};
        std::vector<std::function<T (const std::vector<T>& x, const T t)>> _ode_system;


        Test1(){
            T
                k = 20,
                m = 0.3,
                w = k/m;
            
            // dx/dt
            auto func_1 = [](const std::vector<T>& x, const T t) -> T{
                return x[1];
            };


            // dy/dt
            auto func_2 = [w](const std::vector<T>& x, const T t) -> T{
                return - w * x[0];
            };

            _ode_system.push_back(func_1);
            _ode_system.push_back(func_2);
        }

};
