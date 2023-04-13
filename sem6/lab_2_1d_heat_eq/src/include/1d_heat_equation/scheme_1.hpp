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
    fout << ",u" << x.size() - 1;

    while (x[x.size() - 1] < end_x)
    {   
        x.push_back(x[x.size() - 1] + h);
        fout << ",u" << x.size() - 1;

        // if (x[x.size() - 1] + h > end_x){
        //     break;
        // }

    }
    
    fout << "\n" << start_time;
    for (std::size_t i = 0; i < x.size(); ++i) {
        fout << "," << x[i];
    }
    fout << "\n" << start_time;
    y.push_back(func[1](x[0], start_time));
    fout << "," << y[0];
    for (std::size_t i = 1; i < x.size()-1; ++i) {
        y.push_back(func[0](x[i], start_time));
        fout << "," << y[i];
    }
    y.push_back(func[2](x[x.size()-1], start_time));
    fout << "," << y[x.size()-1];

    std::vector<T> tmp(y);
    std::size_t counter = 0;

    T C = 0;
    C = tau / (h * h);
    T ai, ai1;
    while (start_time + tau <= end_time){
        start_time += tau;
        ++counter;
        if(counter > y.size() / 2){
            break;
        }
        for(std::size_t i = counter; i < tmp.size() - counter; ++i){
            ai1 = func[3](tmp[i+1] - h / 2, start_time);
            ai = func[3](tmp[i] - h / 2, start_time);
            y[i] = C * ( ai1 * (tmp[i + 1] - tmp[i]) - ai * (tmp[i] - tmp[i - 1])) + tmp[i];
        }
        
        // for(std::size_t i = 1; i < y.size()-1; ++i){
        //     ai1 = func[3](tmp[i+1] - h / 2, start_time);
        //     ai = func[3](tmp[i] - h / 2, start_time);
        //     y[i] = C * ( ai1 * (tmp[i + 1] - tmp[i]) - ai * (tmp[i] - tmp[i - 1])) + tmp[i];
        // }

        y[0] = func[1](x[0], start_time);
        y[y.size()-1] = func[2](x[x.size()-1], start_time);

        fout << "\n" << start_time;
        for (std::size_t i = 0; i < y.size(); ++i) {
            fout << "," << y[i];
        }
        tmp.assign(y.begin(), y.end());
    }
};
