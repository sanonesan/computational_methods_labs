#pragma once

#include<vector>
#include<cmath>
#include<functional> 

/*
// ...........Тест..K = const.............. //
*/
template<class T>
class Test1{

    public:

        // std::vector<T> _x0 = {1.,0};
        std::vector<std::function<T (const T x, const T t)>> _system;

        Test1(){

            T 
                c = 500., 
                p = 7800.;
            
            auto K = [](const T x, const T t) -> T{
                return 45.4;
            };

            // dx/dt
            auto u0 = [](const T x, const T t) -> T{
                return - exp(sin(x) * sin(x)) * cos(t);
            };

            _system.push_back(u0);
            _system.push_back(K);
            
        }

};