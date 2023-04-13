#pragma once

#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <vector>

template <typename T, typename F>
void scheme1(T start_time, T end_time, T tau, T start_x, T end_x, T h, const std::vector<F> &func, const std::string &out_path) {

    std::ofstream fout(out_path);
    if (!fout) {
        std::cout << "\n error \n";
        return;
    }

    fout << std::scientific;
	fout << std::setprecision(8);

    fout << "time";

    std::vector<T> x;
    std::vector<T> y;
    x.push_back(start_x);
    while (x[x.size() - 1] < end_x)
    {   
        fout << ",u" << x.size() - 1;
        if (x[x.size() - 1] + h > end_x){
            break;
        }
        x.push_back(x[x.size() - 1] + h);

    }
    
    fout << "\n" << start_time;
    for (std::size_t i = 0; i < x.size(); ++i) {
        fout << "," << x[i];
    }
    fout << "\n" << start_time;
    for (std::size_t i = 0; i < x.size(); ++i) {
        y.push_back(func[0](x[i], start_time));
        fout << "," << y[i];
    }

    std::vector<T> tmp(y);
    std::size_t counter;
    while (start_time <= end_time){
        start_time += tau;
        ++counter;
        if(counter > y.size() / 2){
            break;
        }
        for(std::size_t i = counter; i < tmp.size() - counter; ++i){
            y[i] = tau / h / h * (tmp[i + 1] - 2 * tmp[i] + tmp[i - 1] ) + tmp[i];
        }
        // for(std::size_t i = 0; i < x.size(); ++i){
        //     y[i] = tau / h / h * (tmp[i + 1] - 2 * tmp[i] + tmp[i - 1]) + tmp[i];
        // }
        fout << "\n" << start_time;
        for (std::size_t i = 0; i < x.size(); ++i) {
            fout << "," << y[i];
        }
        tmp.assign(y.begin(), y.end());
    }
};
