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
        T x0 = 0.;
        T xl = 1.;
        T u0 = 0.;

        Test1(){

            T 
                c = 500., 
                p = 7800.;
            
            auto K = [](T x, const T t) -> T{
                return 1.;
            };

            auto u_x_0 = [this](const T x, const T t) -> T{
                return this->u0 + x * (this->xl - x);
            };

            auto u_0_t = [this](const T x, const T t) -> T{
                return this->u0;
            };

            auto u_L_t = [this](const T x, const T t) -> T{
                return this->u0;
            };

            _system.push_back(u_x_0);
            _system.push_back(u_0_t);
            _system.push_back(u_L_t);
            _system.push_back(K);
            
        }

};