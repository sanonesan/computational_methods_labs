#include<iostream>
#include<fstream>
#include<vector>
#define eps 1e-9
#include "unlinear_methods_test_functions.cpp"

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;

ofstream file1;
ofstream file2;
ofstream file3;
ofstream file4;
ofstream file5;
ofstream file6;
ofstream file7;
ofstream file8;
ofstream file9;
ofstream file10;
string path = "../res/";

//Вывод вектора
template<typename T>
void print_vec(const vector<T>& vec)
{
    cout << "( ";
	for (size_t i = 0; i < vec.size() - 1; i++)
		cout << vec[i] << "\t";

	cout << vec[vec.size() - 1] ;
	cout << " )^T\n";
}

//Вывод вектора пар локализации
template<typename T>
void print_vec(const vector<pair<T, T>>& vec)
{   cout << "( ";
	for (size_t i = 0; i < vec.size() - 1; i++)
		cout << "{ " << vec[i].first << "\t" << vec[i].second << " } \t";

	cout << "{ " << vec[vec.size() - 1].first << "\t" << vec[vec.size() - 1].second << " }";    
	cout << " ) \n";
}


// умножение матрицы и столбца
template<typename T>
vector<T> mult(const matrix<T>& A1, const vector<T>& A2) {
	size_t n = A2.size();
	vector<T> res;
	res.reserve(n); // оптимизируем память
	for (size_t i = 0; i < n; ++i) {
		T c = 0;
		for (size_t j = 0; j < n; ++j) {
			c += A1[i][j] * A2[j];
		}
		res.push_back(c);

	}
	return res;
}

// норма (евклидова) вектора
template<typename T>
T norm(const vector<T>& x) {
	T res = 0;
	size_t n = x.size();
	for (size_t i = 0; i < n; ++i)
		res += x[i] * x[i];
	return sqrt(res);
}


// Обращение матрицы 2x2 (просто ввел формулы из Wolfram Mathematica)
template<typename T>
matrix<T> inv(const matrix<T>& m) {
	T det = m[0][0] * m[1][1] - m[1][0] * m[0][1];
	return matrix<T> {
		{  m[1][1] / det, -m[0][1] / det },
		{ -m[1][0] / det,  m[0][0] / det }
	};
}


//error order
template<typename T>
void error_order(vector<T> &x_mas, T x){
    for (size_t i = 0; i < x_mas.size() - 2; ++i){
		cout << "Порядок p = " << log(abs((x_mas[i+2] - x) / (x_mas[i+1] - x))) / log(abs((x_mas[i+1] - x) / (x_mas[i] - x))) << "\n";
    }
}

// Численная производная функции
template<typename T, typename F>
auto num_deriv(F& foo, T _eps) {
	return [&foo, _eps](T x) {
		return (foo(x + _eps) - foo(x)) / _eps;
	};
}

// Вывод отрезков локализации
template<typename T>
ostream& operator<<(ostream& out, vector<pair<T, T>> intervals) {
	if (intervals.size() > 0) {
		pair<T, T> inter = intervals[0];
		out << '[' << inter.first << ", " << inter.second << "]";
	}
	for (size_t i = 1; i < intervals.size(); ++i) {
		pair<T, T> inter = intervals[i];
		out << "; [" << inter.first << ", " << inter.second << "]";
	}
	return out;
}


//////////////////////////////////////////////////
////* unlinear_eq_methods METHODS REALIZATION*////
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
        //cout << grid[i] << "\n";
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
    error_order(x_mas, xk);
    x_mas.clear();

    return (fabs(xk) > eps) ? xk : 0.0;    
}


template<typename T, typename Func>
T  chord_method(Func foo,  pair<T,T> section, size_t& iter_counter){
    
    T a = section.first;
	T b = section.second;
	T xk;

	vector<T> x_mas;
	x_mas.push_back(a);
	x_mas.push_back(b);

    while(true){

        xk = a - (b - a) * foo(a) / (foo(b) - foo(a));
        b = a;
        a = xk;
        x_mas.push_back(xk);

        if (fabs(a - b) < eps){
            break;
        }
    }
    xk = a;

	error_order(x_mas, xk);
    x_mas.clear();
	return xk;
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
        pair<T,T> temp;

        xk = x - foo(x) / diff_foo(x);
        
        if (xk > section.second){
            temp = {section.first, xk};
            xk = chord(foo, temp);
        }
        if (xk < section.first){
            temp = {xk, section.second};
            xk = chord(foo, temp);
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

	size_t iter_counter = 0;

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

		++iter_counter;

	} while (coef > eps);

	cout << "Количество итераций: " << iter_counter << "\n";

	return xk;
}

template<typename T, typename Func, typename J>
vector<T>  newton_sys_method1(Func foo1, Func foo2, J jacobi, const vector<T>& x0, size_t &iter_counter)
{
	vector<T> xk = x0;
	T coef;

	iter_counter = 0;

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

		++iter_counter;

	} while (coef > eps);

	//cout << "Количество итераций: " << iter_counter << "\n";

	return xk;
}



template<typename T, typename Func, typename J>
vector<T>  newton_sys_method_mod(Func foo1, Func foo2, J jacobi, pair<T, T> section1, pair<T, T> section2, const vector<T>& x0)
{
	vector<T> xk = x0;
	T coef;

    pair<T,T> temp;

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
        
        T f_bx = foo1(section1.second, newXk[1]);
        T f_ax = foo1(section1.first, newXk[1]);
        T f_xk = foo1(newXk[0], newXk[1]);

        // return (foo(section.first) * section.second - foo(section.second) * section.first)
		// / (foo(section.first) - foo(section.second));   
		if (newXk[0] < section1.first){
            newXk[0] = (f_ax * newXk[0] - f_xk * section1.first ) / (f_ax - f_xk);
        }
        if (newXk[0] > section1.second){
            newXk[0] =(f_xk * newXk[0] - f_bx * section1.second ) / (f_xk - f_bx);
        }

        T f_ay = foo2(newXk[0], section2.first);
        T f_by = foo2(newXk[0], section2.second);
        T f_yk = foo2(newXk[0], newXk[1]);

        if (newXk[1] < section2.first){
            newXk[1] = (f_ay * newXk[1] - f_xk * section2.first ) / (f_ay - f_yk);
        }
        if (newXk[1] > section2.second){
            newXk[1] = (f_yk * newXk[1] - f_by * section2.second ) / (f_yk - f_by);
        }        

		coef = norm(jacobiInvMultF);

		xk = newXk;

	} while (coef > eps);

	cout << "Количество итераций: " << iter_counter << "\n";

	return xk;
}

//////////////////////////////////////////////////
////* unlinear_eq_methods METHODS REALIZATION*////
//////////////////////////////////////////////////


void foo(size_t k, double x, double y, string path);

int main(int args, char **argv){

    typedef double Type;

    vector<Type> x0test4{ -2.433, -9.654 };
	vector<Type> x0test5{ -1., 4. };
	vector<Type> x0test6{ 0.9, 0.5 };

    size_t iter = 0;
    pair<Type, Type> a = {0, 1};
    double x = 0;

    cout << "\n/--------------------------------------------/ \n";
    cout << "/--------------------------------------------/ \n";

    cout << "\nTEST0: \n";

    auto testing = test0<Type>;
    
    cout << "Sections with roots: \n";
    vector<pair<Type, Type>> intervals = localization(testing, a, 100);
    print_vec(intervals);

    cout << "\nBisection: \n";
    vector<Type> X;
    for(size_t i = 0; i < intervals.size(); ++i){
        iter = 0;
        X.push_back(bisection(testing, intervals[i], iter));
    }
    print_vec(X);
    X.clear();

    cout << "\nChords: \n";
    for(size_t i = 0; i < intervals.size(); ++i){
        iter = 0;
        X.push_back(chord_method(testing, intervals[i], iter));
    }
    print_vec(X);
    X.clear();


    cout << "\nNewton_method: \n";
    for(size_t i = 0; i < intervals.size(); ++i){
        iter = 0;
        X.push_back(newton_method(testing, test0_diff<Type>, intervals[i], iter));
    }
    print_vec(X);

    vector<Type> Y;
    for(size_t i = 0; i < intervals.size(); ++i){
        iter = 0;
        Y.push_back(newton_method(testing, num_deriv(test0<Type>, 1e-3), intervals[i], iter));
    }
    print_vec(Y);

    X.clear();
    Y.clear();

    cout << "\n/--------------------------------------------/ \n";
    cout << "/-----------------СИСТЕМЫ--------------------/ \n";
    cout << "/--------------------------------------------/ \n";

	cout << "\nКлассический метод Ньютона (система) ТЕСТ 4\n";
	vector<Type> rootNewSys = newton_sys_method(
        test4f1<Type>, 
        test4f2<Type>, 
        test4Jacobi<Type>,
        //vector<Type>{1.5, 4.3}
		x0test4
        );
	print_vec(rootNewSys);

    cout << "\nМодифицированный метод Ньютона (система) ТЕСТ 4\n";
	rootNewSys = newton_sys_method_mod(
        test4f1<Type>, 
        test4f2<Type>, 
        test4Jacobi<Type>,
        pair<Type, Type>{-10, 10},
        pair<Type, Type>{-10, 10},
		x0test4
        );
	print_vec(rootNewSys);


    cout << "\nКлассический метод Ньютона (система) ТЕСТ 5\n";
    rootNewSys = newton_sys_method(
        test5f1<Type>, 
        test5f2<Type>, 
        test5Jacobi<Type>,
		x0test5
        );
	print_vec(rootNewSys);

	cout << "\nМодифицированный метод Ньютона (система) ТЕСТ 5\n";
    rootNewSys = newton_sys_method_mod(
        test5f1<Type>, 
        test5f2<Type>, 
        test5Jacobi<Type>,
        pair<Type, Type>{-10, 10},
        pair<Type, Type>{-10, 10},
		x0test5
        );
	print_vec(rootNewSys);

    
    file1.open(path + "1_4.txt");
	file2.open(path + "5_6.txt");
	file3.open(path + "7_8.txt");
	file4.open(path + "9_10.txt");
	file5.open(path + "11_12.txt");
	file6.open(path + "13_15.txt");
	file7.open(path + "16_20.txt");
	file8.open(path + "21_25.txt");
	file9.open(path + "26_30.txt");
	file10.open(path + "30.txt");


    Type i = -10, j = -10;
	Type h = 0.2;
	vector<Type> test;
    iter = 0;
	while (i < 10)
	{
		while (j < 10) {
			rootNewSys = newton_sys_method1(test4f1<Type>, test4f2<Type>, test4Jacobi<Type>, vector<Type>{i, j}, iter);
            //cout << i << "\t" << j << "\t";
			foo(iter, i, j, path);
			j += h;
		}
		j = -10;
		i += h;
	}


    cout << "\n/--------------------------------------------/ \n";
    cout << "/------------------ОТЧЕТ---------------------/ \n";
    cout << "/--------------------------------------------/ \n\n";

	cout << fixed;
	cout.precision(15);
	auto f6 = test6<Type>;
	auto f6_diff = test6_diff<Type>;
	auto f6_diff_num = num_deriv(f6, eps);
    double exact_root = 0.608246;
    pair<Type, Type> section{0., 1.};

    cout << "Метод бисекции:\n";
	iter = 0;
	Type root = bisection(f6, section, iter);
	cout << "Количество итераций: " << iter << "\n";
	cout << "Результат: " << root << "\n";
	cout << "Достигнутая точность: " << abs(exact_root - root) << "\n";
	cout << "Невязка: " << abs(f6(root)) << "\n";

    cout << "\n/--------------------------------------------/ \n\n";


    cout << "Метод хорд:\n";
	iter = 0;
	root = chord_method(f6, section, iter);
	cout << "Количество итераций: " << iter << "\n";
	cout << "Результат: " << root << "\n";
	cout << "Достигнутая точность: " << abs(exact_root - root) << "\n";
	cout << "Невязка: " << abs(f6(root)) << "\n";

    cout << "\n/--------------------------------------------/ \n\n";

    cout << "Метод Ньютона (аналитическая производная):\n";
	iter = 0;
    root = newton_method(f6, f6_diff, section, iter);
	cout << "Количество итераций: " << iter << "\n";
	cout << "Результат: " << root << "\n";
	cout << "Достигнутая точность: " << abs(exact_root - root) << "\n";
	cout << "Невязка: " << abs(f6(root)) << "\n";

    cout << "\n/--------------------------------------------/ \n\n";

    cout << "Метод Ньютона (численная производная):\n";
	iter = 0;
    root = newton_method(f6, f6_diff_num, section, iter);
	cout << "Количество итераций: " << iter << "\n";
	cout << "Результат: " << root << "\n";
	cout << "Достигнутая точность: " << abs(exact_root - root) << "\n";
	cout << "Невязка: " << abs(f6(root)) << "\n";
    


    return 0;
}




void foo(size_t k, double x, double y, string path)
{
    
	if (k < 5) {
		file1 << x << "  " << y << "\n"; 
        return;
	}
	if (k < 7) {
		file2 << x << "  " << y << std::endl; return;
	}
	if (k < 9) {
		file3 << x << "  " << y << std::endl; return;
	}
	if (k < 11) {
		file4 << x << "  " << y << std::endl; return;
	}
	if (k < 13) {
		file5 << x << "  " << y << std::endl; return;
	}
	if (k < 16) {
		file6 << x << "  " << y << std::endl; return;
	}
	if (k < 21) {
		file7 << x << "  " << y << std::endl; return;
	}
	if (k < 26) {
		file8 << x << "  " << y << std::endl; return;
	}
	if (k < 31) {
		file9 << x << "  " << y << std::endl; return;
	}
	if (k > 30) {
		file10 << x << "  " << y << std::endl; return;
	}
	
}