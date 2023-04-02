#pragma once

#include<vector>
#include<cmath>
#include<functional> 

/*
// ...........Тест..из..методички.............. //
*/
template<typename T>
class Test02{

    public:

        std::vector<T> _x0 = {1.};
        std::vector<std::function<T (const std::vector<T>& x, const T t)>> _ode_system;


        Test02(){
            
            // dx/dt
            auto func_1 = [](const std::vector<T>& x, const T t) -> T{
                return 2 * exp(sin(t) * sin(t)) * cos(t) * sin(t);
            };

            _ode_system.push_back(func_1);
            // _ode_system.push_back(func_2);
        }

};
