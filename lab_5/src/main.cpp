#include<iostream>
#include<fstream>
#include<vector>
#define eps 1e-9
//#include "unlinear_eq_methods.h"
#include "unlinear_methods_test_functions.cpp"
#include "helpfunc.h"

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;



//////////////////////////////////////////////////
/* CLASS unlinear_eq_methods METHODS REALIZATION*/
//////////////////////////////////////////////////


template<typename T, typename Func>
vector<pair<T,T>>  localization(Func foo, pair<T,T> main_section, size_t n){

    T h = (main_section.second - main_section.first) / (n-1);
    vector<T> grid;
    for(size_t i = 0; i < n; ++i){
        grid.push_back(main_section.first + h * i);
    }

    pair<T,T> section_with_root;
    vector<pair<T,T>> localized_roots;

    for(size_t i = 0; i < n; ++i){
        //cout << grid[i] << endl;
        if(foo(grid[i]) * foo(grid[i+1]) < 0.0){
            section_with_root.first = grid[i];
            section_with_root.second = grid[i+1];
            localized_roots.push_back(section_with_root);
        }
    }
    // if(localized_roots.size() == 0){
    //     throw("no roots in this section");
    // }
    return localized_roots;
}


template<typename T, typename Func>
T  bisection(Func foo, pair<T,T> section, size_t& iter_counter){

    T ak = section.first;
    T bk = section.second;
    T f_ak = foo(ak);
    T f_bk = foo(bk);
    T xk = (ak + bk) / 2;
    T f_xk = foo(xk);
    vector<T> x_mas;

    while(true){

        ++iter_counter;

        if (f_ak * f_xk <= 0){
            bk = xk;
            f_bk = f_xk;
        }
        else{
            ak = xk;
            f_ak = f_xk;
        }        

        xk = (ak + bk) / 2;
        f_xk = foo(xk);
        
        x_mas.push_back(xk);

        if(fabs(ak - bk) / 2 < eps)
            break;
        
    }
    //error_order(x_mas, 1.0);
    x_mas.clear();

    return xk;    
}


template<typename T, typename Func>
T  chord_method(Func foo,  pair<T,T> section, size_t& iter_counter){
    
    T xk = section.first;
	T xk1 = section.second;
	T xk2;

	vector<T> x_mas;
	x_mas.push_back(xk);
	x_mas.push_back(xk1);

	while(true){

        ++iter_counter;
		xk2 = xk - f(xk) * (xk1 - xk) / (f(xk1) - f(xk));
		x_mas.push_back(xk2);
		xk = xk1;
		xk1 = xk2;

		if (fabs(xk2 - xk1) < eps)
            break;

    }
	//error_order(x_mas, 1.0);
    x_mas.clear();
	return xk2;
}

template<typename T, typename Func>
T  chord(Func foo,  pair<T,T> section){
    
    return (foo(section.first) * section.second - foo(section.second) * section.first)
		/ (foo(section.first) - foo(section.second));   

}


template<typename T, typename Func, typename Diff>
T  newton_method(Func foo, Diff diff_foo, pair<T,T> section, size_t& iter_counter){

    T x, xk = chord(foo, section);
    vector<T> temp;

    do {
        ++iter_counter;
        
        x = xk;
        temp.push_back(x);
        
        xk = x - foo(x) / diff_foo(x);
        
        if (xk > section.second){
            xk = chord(foo, {section.first, xk});
        }
        if (xk < section.first){
            xk = chord(foo, {xk, section.second});
        }

    } while (fabs(xk - x) > eps);

    return x;    
}



// Классический метод Ньютона для системы уравенений
// На каждой итерации считается матрица Якоби
template<typename T, typename Func, typename J>
vector<T>  newton_sys_method(Func foo1, Func foo2, J jacobi, const vector<T>& x0)
{
	vector<T> xk = x0;
	T coef;

	size_t iterCount = 0;

	do {

		matrix<T> jacobiInv = inv(jacobi(xk[0], xk[1]));

		vector<T> F{ foo1(xk[0], xk[1]), foo2(xk[0], xk[1]) };

		vector<T> jacobiInvMultF = mult(jacobiInv, F);

		vector<T> newXk;

		newXk.reserve(x0.size());
		for (size_t i = 0; i < x0.size(); ++i) {
			newXk.push_back(xk[i] - jacobiInvMultF[i]);
		}

		coef = norm(jacobiInvMultF);

		xk = newXk;

		iterCount++;

	} while (coef > eps);

	cout << "Количество итераций: " << iterCount << endl;

	return xk;
}



template<typename T, typename Func, typename J>
vector<T>  newton_sys_method_mod(Func foo1, Func foo2, J jacobi, pair<T, T> section1, pair<T, T> section2, const vector<T>& x0)
{
	vector<T> xk = x0;
	T coef;

	size_t iter_counter = 0;

	do {

		++iter_counter;
		// F'^-1
		matrix<T> jacobiInv = inv(jacobi(xk[0], xk[1]));
		vector<T> F{ foo1(xk[0], xk[1]), foo2(xk[0], xk[1]) };

		// F'^-1 * F
		vector<T> jacobiInvMultF = mult(jacobiInv, F);

		vector<T> newXk = xk;
		newXk[0] -= jacobiInvMultF[0];
		newXk[1] -= jacobiInvMultF[1];

        // Для защиты от выхода за область применяем гибридный метод:
		// Внешние итерации - одна итерация нелинейного метода Зейделя
		// Внутренние - одна итерация метода хорд

		if (newXk[0] < section1.first){
            newXk[0] = chord(foo1, {section1.first, newXk[0]});
        }
        if (newXk[0] > section1.second){
            newXk[0] = chord(foo1, {newXk[0], section1.second});
        }

        if (newXk[1] < section2.first){
            newXk[1] = chord(foo2, {section2.first, newXk[1]});
        }
        if (newXk[1] > section2.second){
            newXk[1] = chord(foo2, {newXk[1], section2.second});
        }        

		coef = norm(jacobiInvMultF);

		xk = newXk;

	} while (coef > eps);

	cout << "Количество итераций: " << iter_counter << endl;

	return xk;
}

//////////////////////////////////////////////////
/* CLASS unlinear_eq_methods METHODS REALIZATION*/
//////////////////////////////////////////////////
template<typename T>
void print_vec(const vector<pair<T, T>>& vec)
{
	for (size_t i = 0; i < vec.size(); i++)
		cout << vec[i].first << "\t" << vec[i].second << endl;
	cout << endl;
}


int main(int args, char **argv){

    typedef double Type;

    //unlinear_eq_methods solver;

    size_t iter = 0;
    pair<Type, Type> a = {-10.0, 10};
    
   // double x = bisection(test1<Type>, a, iter);
    
    //x = solver.bisection(test1<Type>, a, iter);
    // size_t k = 15;
    auto intervals = localization(test0<Type>, a, 100);

    cout << intervals.size() << endl;
    print_vec(intervals);
    //print_vec(intervals);


    return 0;
}