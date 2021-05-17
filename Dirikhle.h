#pragma once
// Dirikhle_system.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <functional>


const double eps = std::numeric_limits<double>::epsilon();


double testFunc(double x) {
    return x <= 0 ? x * x * x + 3 * x * x : -x * x * x + 3 * x * x;
}
double testFuncDiv1(double x) {
    return x <= 0 ? 3 * x * x + 6  * x : -3* x * x + 6  * x;
}
double testFuncDiv2(double x) {
    return x <= 0 ? 6  * x + 6  : -6 * x + 6 ;
}
double mainFunc(double x) {
    return sin(cos(x));
}
double mainFuncDiv1(double x) {
    return -sin(x) * cos(cos(x));
}
double mainFuncDiv2(double x) {
    return -sin(x)*sin(x)*sin(cos(x)) - cos(x)*cos(cos(x));
}

double mainFunc1(double x) {
    return sin(cos(x)) + cos(10*x);
}
double mainFunc1Div1(double x) {
    return -sin(x) * cos(cos(x)) - 10*sin(10*x);
}
double mainFunc1Div2(double x) {
    return -sin(x)*sin(x)*sin(cos(x)) - cos(x)*cos(cos(x)) - 100*cos(10*x);
}

struct TResults {
    std::vector<std::pair<double, double> > res_vec;
    std::vector<double> local_mistake_vec;
    std::vector<double> h_vec;
    std::vector<double> a;
    std::vector<double> d;
    std::vector<double> b;
    std::vector<double> c;
    uint64_t ND = 0;
    uint64_t NH = 0;
    //double max_local_mistake;
};

void straightRun(const  std::vector <double>& longDiag, const  std::vector <double>& shortDiagUp, const  std::vector <double>& shortDiagDown, const  std::vector<double>& answer, std::vector<double>& alphas, std::vector <double>& betas) {
    alphas.push_back(-shortDiagUp.at(0) / longDiag.at(0));
    betas.push_back(answer.at(0) / longDiag.at(0));
    for (size_t i = 1; i < longDiag.size() - 1; i++) {
        alphas.push_back(-shortDiagUp.at(i) / (longDiag.at(i) + shortDiagDown.at(i - 1) * alphas.at(i - 1)));
        betas.push_back((answer.at(i) - shortDiagDown.at(i - 1) * betas.at(i - 1)) / (longDiag.at(i) + shortDiagDown.at(i - 1) * alphas.at(i - 1)));
    }
    betas.push_back((answer.back() - shortDiagDown.back() * betas.back()) / (longDiag.back() + shortDiagDown.back() * alphas.back()));
    alphas.push_back(double(0));
}
void backRun(double an, double cn, double answern, const  std::vector <double>& alphas, const  std::vector <double>& betas, TResults& res) {
    res.res_vec[res.res_vec.size() - 1].second = betas.back();
    size_t sz = alphas.size();
    for (size_t i = 1; i < sz; i++) {
        res.res_vec[res.res_vec.size() - 1 - i].second = alphas.at(sz - i - 1) * res.res_vec[res.res_vec.size() - i].second + betas.at(sz - i - 1);
    }
}
void  tridiagonalMatrixAlg(const  std::vector <double>& longDiag, const  std::vector <double>& shortDiagUp, const  std::vector <double>& shortDiagDown, const  std::vector <double>& answer, TResults& res) {
    if (shortDiagUp.size() != shortDiagDown.size() || longDiag.size() == shortDiagDown.size()) {
        throw 1;
    }
    else {
        std::vector <double> alphas;
        std::vector <double> betas;
        straightRun(longDiag, shortDiagUp, shortDiagDown, answer, alphas, betas); // Прямой ход 
        std::vector<double> x;
        backRun(shortDiagDown.back(), longDiag.back(), answer.back(), alphas, betas, res); // Обратный ход
    }

}

TResults tfunc1(double start, double b, int n, std::function<double(double)> f) {
    TResults result;
    double h = (b-start) / n;

    for (int i = 0; i < n+1; i++) {
        result.res_vec.push_back({ 0, 0 });
    }
    result.a.reserve(n);
    result.b.reserve(n);
    result.c.reserve(n + 1);
    result.d.reserve(n);
    std::vector<double> longDiag;
    std::vector<double> shortUP;
    std::vector<double> shortDown;
    std::vector<double> answer;

    answer.push_back(0);
    longDiag.push_back(1);
    shortUP.push_back(0);

    for (int i = 0; i < n-1; i++) {
        shortDown.push_back(h);
        longDiag.push_back(4 * h);
        shortUP.push_back(h);
        answer.push_back(6 * (((f(start + (i + 2) * h) - f(start + (i + 1) * h)) / h) - ((f(start + (i + 1) * h) - f(start + i * h)) / h)));
    }

    shortDown.push_back(0);
    longDiag.push_back(1);
    answer.push_back(0);
    tridiagonalMatrixAlg(longDiag, shortUP, shortDown, answer, result);

    result.res_vec[0].first = start;
    for (int i = 0; i < n; i++) {
        result.a.push_back( f(start + (i + 1) * h));
        result.res_vec[i + 1].first = start + (i + 1) * h;
        result.c.push_back( result.res_vec[i].second);
        result.d.push_back( (result.res_vec[i + 1].second - result.res_vec[i].second)/h);
        result.b.push_back ( (f(start + (i + 1) * h) - f(start + (i)*h)) / h + result.res_vec[i + 1].second * h / 3 + result.res_vec[i].second * h / 6);
    }
    result.c.push_back(result.res_vec[n].second);

    return result;
}

std::string doubleToString(double num) {

    std::stringstream ss;
    ss << num;

    return ss.str();
}

double setPresision(double num, int presision) {
    num *= pow(10, presision + 1);

    return double(int(num)) / pow(10, presision + 1);
}


std::string writeFinalTable(int n, int m, std::vector<std::vector<double>> v, double a, double b, double c, double d, int S, int S_2, double eps, double eps_max, double eps_max_2, double error_max, double r_norm, double w, double accuracy) {
    std::ofstream outfile("out.txt", std::ofstream::app);
    outfile << "Таблица № " << 1 << std::endl;
    for (size_t j = 0; j < m + 3; j++) {
        if (j == 0) {
            outfile << std::setw(5) << "|" << std::setw(7) << "i" << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << i << "|";
            }

            outfile << std::endl;
        } else if (j == 1) {
            outfile << std::setw(4) << "j" << "|" << std::setw(7) << "Y/X" << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << std::setprecision(5) << a + i * (b - a) / n << "|";
            }
            outfile << std::endl;
        } else {
            outfile << std::setw(4) << j - 2 << "|" << std::setw(7) << std::setprecision(4) << c + (j - 2) * (d - c) / m << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << std::setprecision(5) << v[i][j - 2] << "|";
            }
            outfile << std::endl;
        }
    }
    outfile << "\nТаблица № " << 7 << "\n\n" << std::endl;
    std::string returned_string = std::string("Параметр W: ") + doubleToString(w) + "\n\nЧисло шагов: " + doubleToString(S) + "\n\nЧисло шагов с двойным шагом: " + doubleToString(S_2) + "\n\nТочность метода: " +
        doubleToString(eps) + "\n\nДостигнутая точность: " + doubleToString(eps_max) + "\n\nДостигнутая точность с двойным шагом: " + doubleToString(eps_max_2) + "\n\nНорма невязки: " + doubleToString(r_norm) + "\n\nТочность max |v - v_2|: " + doubleToString(accuracy);
    
    //outfile << std::endl << "параметр w: " << w << std::endl;
    //outfile << std::endl << "число шагов: " << S << std::endl;
    //outfile << std::endl << "число шагов с двойным шагом: " << S_2 << std::endl;
    //outfile << std::endl << "точность метода: " << eps << std::endl;
    //outfile << std::endl << "достигнутая точность: " << eps_max << std::endl;
    //outfile << std::endl << "достигнутая точность с двойным шагом: " << eps_max_2 << std::endl;
    ////outfile << std::endl << "погрешность: " << error_max << std::endl;
    //outfile << std::endl << "норма невязки: " << r_norm << std::endl;
    //outfile << std::endl << "точность max|v-v_2|: " << accuracy << std::endl;
    outfile << returned_string;
    outfile.close();
    return returned_string;
}

void writeTable(int n, int m, std::vector<std::vector<double>> v, double a, double b, double c, double d, int tableNumber) {
    std::ofstream outfile("out.txt", std::ofstream::app);
    outfile << "Таблица № " << tableNumber << std::endl;
    for (size_t j = 0; j < m + 3; j++) {
        if (j == 0) {
            outfile << std::setw(5) << "|" << std::setw(7) << "i" << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << i << "|";
            }

            outfile << std::endl;
        } else if (j == 1) {
            outfile << std::setw(4) << "j" << "|" << std::setw(7) << "Y/X" << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << std::setprecision(5) << a + i * (b - a) / n << "|";
            }
            outfile << std::endl;
        } else {
            outfile << std::setw(4) << j - 2 << "|" << std::setw(7) << std::setprecision(4) << c + (j - 2) * (d - c) / m << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << std::setprecision(5) << v[i][j - 2] << "|";
            }
            outfile << std::endl;
        }
    }
    outfile << std::endl;
    outfile.close();
}

void writeHeader(int n, int m) {
    std::ofstream outfile("out.txt");

    outfile << "\n\nПрименение итерационного метода для решения разностных схем\nна примере задачи Дирихле для уравнения Пуассона\nКоманда №7." << std::endl;
    outfile << "\nu = exp(sin(pi * x * y) * sin(pi * x * y)). Сетка " << n << "x" << m << ". Границы: a=0, b=2, c=0, d=1\n" << std::endl;
}

template <typename Function>
double calculateAlpha(double h2,double k2,double a, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& hs, Function func,const int n, const int m,double h,double k,double ag,double cg) {

    double top = 0;
    double value=0;
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            value= (a * v[i][j] + h2*v[i-1][j]+h2*v[i+1][j]+k2*v[i][j-1]+k2*v[i][j+1] - func(ag+i*h,cg+j*k))*hs[i][j];
            top += value;
        }
    }

    double bottom = 0;
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }

    return (-top / bottom);
}


template <typename Function>
double calculateBeta(double h2, double k2, double a, std::vector<std::vector<double>>& r, std::vector<std::vector<double>>& hs, Function func, const int n, const int m, double h, double k, double ag, double cg) {

    double top = 0;
    double value;
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1] ) * r[i][j];
            top += value;
        }
    }

    double bottom = 0;
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }

    return top / bottom;
}

template <typename Function>
double calculateAlphaCut(double h2, double k2, double a, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& hs, Function func, const int n, const int m, double h, double k, double ag, double cg) {

    double top = 0;
    double value = 0;
    double bottom = 0;
   
    // 1
    for (int i = 1; i <= n / 4; i++) {
        for (int j = 1; j < 3 * m / 4; j++) {
       
            value = (a * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(ag + i * h, cg + j * k)) * hs[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }
    //2
    for (int i = n / 4 + 1; i < n / 2; i++) {
        for (int j = 1; j < m; j++) {
    
            value = (a * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(ag + i * h, cg + j * k)) * hs[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }
    //3
    for (int i = n / 2; i <= 3 * n / 4; i++) {
        for (int j = 3 * m / 4 + 1; j < m; j++) {
     
            value = (a * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(ag + i * h, cg + j * k)) * hs[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }
    //4
    for (int i = 3 * n / 4+1; i < n ; i++) {
        for (int j = 1; j < m; j++) {
    
            value = (a * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(ag + i * h, cg + j * k)) * hs[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }
    //5
    for (int i = n / 2; i <= 3 * n / 4; i++) {
        for (int j = 1; j < m / 4; j++) {
       
            value = (a * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(ag + i * h, cg + j * k)) * hs[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }

    return (-top / bottom);
}


template <typename Function>
double calculateBetaCut(double h2, double k2, double a, std::vector<std::vector<double>>& r, std::vector<std::vector<double>>& hs, Function func, const int n, const int m, double h, double k, double ag, double cg) {

    double top = 0;
    double value;
    double bottom = 0;
  
    // 1
    for (int i = 1; i <= n / 4; i++) {
        for (int j = 1; j < 3 * m / 4; j++) {
      
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * r[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }
    //2
    for (int i = n / 4 + 1; i < n / 2; i++) {
        for (int j = 1; j < m; j++) {
      
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * r[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }
    //3
    for (int i = n / 2; i <= 3 * n / 4; i++) {
        for (int j = 3 * m / 4 + 1; j < m; j++) {
   
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * r[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }
    //4
    for (int i = 3 * n / 4 + 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
         
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * r[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }
    //5
    for (int i = n / 2; i <= 3 * n / 4; i++) {
        for (int j = 1; j < m / 4; j++) {
        
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * r[i][j];
            top += value;
            value = (a * hs[i][j] + h2 * hs[i - 1][j] + h2 * hs[i + 1][j] + k2 * hs[i][j - 1] + k2 * hs[i][j + 1]) * hs[i][j];
            bottom += value;
        }
    }

    return top / bottom;
}


double equalsZero(double num) {
    return (num < eps&& num > -eps) ? 0.0 : num;
}

double w_optimal(double a, double b, double c, double d, double h1, double h2) {
    double l1 = b - a;
    double l2 = d - c;

    double lam_min = (2 * h2 * h2 * sin(M_PI * h1 / (2 * l1)) * sin(M_PI * h1 / (2 * l1)) +
        2 * h1 * h1 * sin(M_PI * h2 / (2 * l2)) * sin(M_PI * h2 / (2 * l2))) / (h1 * h1 + h2 * h2);
    return 2 / (1 + sqrt(lam_min * (2 - lam_min)));
}


template <typename Function>
void solve(std::vector<std::vector<double>>& v, Function func, const int n, const int m, double a, double b, double c, double d, int Nmax, int& S, double& eps, double& eps_max, double& error_max) {
    int i, j; //индексы
    double a2, k2, h2;

    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    std::vector<std::vector<double> >  h(n+1, std::vector<double>(m+1));
    std::vector<std::vector<double> >  accuracy(n + 1, std::vector<double>(m + 1));

    auto m1 = [](double x, double y) {
        return equalsZero(sin(M_PI * y) * sin(M_PI * y));
    };
    auto m2 = [](double x, double y) {
        return equalsZero(sin(M_PI * 2 * y) * sin(M_PI * 2 * y));
    };
    auto m3 = [](double x, double y) {
        return equalsZero(sin(M_PI * x) * sin(M_PI * x));
    };
    auto m4 = [](double x, double y) {
        return equalsZero(sin(M_PI * 2 * x) * sin(M_PI * 2 * x));
    };

    for (j = 0; j < m + 1; j++) {
        v[0][j] = m1(a, c + (d - c) / m * j);
    }
    for (j = 0; j < m + 1; j++) {
        v[n][j] = m2(b, c + (d - c) / m * j);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][0] = m3(a + i * (b - a) / n, c);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][m] = m4(a + i * (b - a) / n, d);
    }



    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            h[i][j] = -(a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b-a)/n, c + j * (d-c)/m));
        }
    }

    double alpha = calculateAlpha(h2, k2, a2, v, h, func, n, m, (b - a) / n, (d - c) / m, a, c);
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < m + 1; j++) {
            v[i][j] += alpha * h[i][j];
        }
    }
  
    eps_max=0;
    double beta = 0;
    S = 1;
    bool flag = false;
   while(!flag) {
        std::vector<std::vector<double> >  r(n+1, std::vector<double>(m+1));
        for (int i = 1; i < n ; i++) {
            for (int j = 1; j < m ; j++) {
              r[i][j] = (a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b-a)/n, c + j * (d-c)/m));
            }
        }
        beta = calculateBeta(h2, k2, a2, r, h, func, n, m, (b - a) / n, (d - c) / m, a, c);

        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m + 1; j++) {
                h[i][j] =   beta * h[i][j] - r[i][j];
            }
        }

        alpha = calculateAlpha(h2, k2, a2, v, h, func, n, m, (b - a) / n, (d - c) / m, a, c);

        double save_component=0;
        eps_max = 0;
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m + 1; j++) {   
                save_component = v[i][j];
                v[i][j] += alpha * h[i][j];
                save_component = abs(v[i][j] - save_component);
                if (save_component > eps_max)
                    eps_max = save_component;

            }
        }
        S += 1;
        if (eps_max <= eps || S >= Nmax) {
            flag = true;
        }

    }
}

template <typename Function>
void solveTest(std::vector<std::vector<double>>& v, Function func, const int n, const int m, double a, double b, double c, double d, int Nmax, int& S, double& eps, double& eps_max, double& error_max) {
    int i, j; //индексы
    double a2, k2, h2;

    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    std::vector<std::vector<double> >  h(n + 1, std::vector<double>(m + 1));
    std::vector<std::vector<double> >  accuracy(n + 1, std::vector<double>(m + 1));

    auto m1 = [](double x, double y) {
        return exp(sin(M_PI * x * y)* sin(M_PI * x * y));
    };
    auto m2 = [](double x, double y) {
        return exp(sin(M_PI * x * y) * sin(M_PI * x * y));
    };
    auto m3 = [](double x, double y) {
        return exp(sin(M_PI * x * y) * sin(M_PI * x * y));
    };
    auto m4 = [](double x, double y) {
        return exp(sin(M_PI * x * y) * sin(M_PI * x * y));
    };

    for (j = 0; j < m + 1; j++) {
        v[0][j] = m1(a, c + (d - c) / m * j);
    }
    for (j = 0; j < m + 1; j++) {
        v[n][j] = m2(b, c + (d - c) / m * j);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][0] = m3(a + i * (b - a) / n, c);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][m] = m4(a + i * (b - a) / n, d);
    }
 

    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j++) {
            h[i][j] = -(a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
        }
    }

    double alpha = calculateAlpha(h2, k2, a2, v, h, func, n, m, (b - a) / n, (d - c) / m, a, c);
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < m + 1; j++) {
            v[i][j] += alpha * h[i][j];
        }
    }

    eps_max = 0;
    double beta = 0;
    S = 1;
    bool flag = false;
    while (!flag) {
        std::vector<std::vector<double> >  r(n + 1, std::vector<double>(m + 1));
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < m; j++) {
                r[i][j] = (a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
            }
        }
        beta = calculateBeta(h2, k2, a2, r, h, func, n, m, (b - a) / n, (d - c) / m, a, c);

        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m + 1; j++) {
                h[i][j] = beta * h[i][j] - r[i][j];
            }
        }

        alpha = calculateAlpha(h2, k2, a2, v, h, func, n, m, (b - a) / n, (d - c) / m, a, c);

        double save_component = 0;
        eps_max = 0;
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m + 1; j++) {
                save_component = v[i][j];
                v[i][j] += alpha * h[i][j];
                save_component = abs(v[i][j] - save_component);
                if (save_component > eps_max)
                    eps_max = save_component;

            }
        }
        S += 1;
        if (eps_max <= eps || S >= Nmax) {
            flag = true;
        }

    }
}


template <typename Function>
void solveCut(std::vector<std::vector<double>>& v, Function func, const int n, const int m, double a, double b, double c, double d, int Nmax, int& S, double& eps, double& eps_max, double& error_max) {
    int i, j; //индексы
    double a2, k2, h2;

    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);

    std::vector<std::vector<double> >  h(n + 1, std::vector<double>(m + 1));
    std::vector<std::vector<double> >  accuracy(n + 1, std::vector<double>(m + 1));

    auto m1 = [](double x, double y) {
        return equalsZero(exp(sin(M_PI*x*y)* sin(M_PI * x * y)));
    };
    

    for (j = 0; j <= 3*m/4 ; j++) {
        v[0][j] = m1(a, c + (d - c) / m * j);
    }
    for (j = 0; j <= m ; j++) {
        v[n][j] = m1(b, c + (d - c) / m * j);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][0] = m1(a + i * (b - a) / n, c);
    }
    for (i = n/4; i < n + 1; i++) {
        v[i][m] = m1(a + i * (b - a) / n, d);
    }

    for (i = n / 2; i <= 3 * n / 4; i++) {
        v[i][m/4] = m1(a+i*(b-a)/n,c+(d-c)/4);
    }

    for (i = n / 2; i <= 3 * n / 4; i++) {
          v[i][3*m/4] = m1(a + i * (b - a) / n, c + 3*(d - c) / 4);
    }

    for (j = m / 4; j <= 3 * m / 4; j++) {
        v[n/2][j] = m1(a+(b-a)/2, c + (d - c) / m * j);
    }

    for (j = m / 4; j <= 3 * m / 4; j++) {
        v[3*n / 4][j] = m1(a+3*(b-a)/4, c + (d - c) / m * j);
    }

    for (i = 1; i <= n/4; i++) {
        v[i][3*m/4] = m1(a + i * (b - a) / n,c+3*(d-c)/4);
    }

    for (j = 3*m / 4; j <=  m ; j++) {
        v[n/4][j] = m1(a + (b - a) / 4, c + (d - c) / m * j);
    }
    std::vector<std::vector<double> >  check(n + 1, std::vector<double>(m + 1));
    // 1
    for (int i = 1; i <= n/4; i++) {
        for (int j = 1; j < 3*m/4; j++) {
            check[i][j] += 1;
            h[i][j] = -(a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
        }
    }
    //2
    for (int i = n/4+1; i < n/2; i++) {
        for (int j = 1; j < m; j++) {
            check[i][j] += 2;
            h[i][j] = -(a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
        }
    }
    //3
    for (int i = n/2; i <= 3*n/4 ; i++) {
        for (int j = 3*m/4+1; j < m; j++) {
            check[i][j] += 3;
            h[i][j] = -(a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
        }
    }
    //4
    for (int i = 3*n/4+1; i < n ; i++) {
        for (int j = 1; j <  m ; j++) {          
            check[i][j] += 4;
            h[i][j] = -(a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
        }
    }
    //5
     for (int i = n/2; i <= 3*n/4; i++) {
        for (int j = 1; j < m/4; j++) { 
            check[i][j] += 5;
            h[i][j] = -(a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
        }
    }

    double alpha = calculateAlphaCut(h2, k2, a2, v, h, func, n, m, (b - a) / n, (d - c) / m, a, c);
   
    // 1
    for (int i = 1; i <= n / 4; i++) {
        for (int j = 1; j < 3 * m / 4; j++) {      
            check[i][j] += 1;
            v[i][j] += alpha * h[i][j];
        }
    }
    //2
    for (int i = n / 4 + 1; i < n / 2; i++) {
        for (int j = 1; j < m; j++) {        
            check[i][j] += 2;
            v[i][j] += alpha * h[i][j];
        }
    }
    //3
    for (int i = n / 2; i <= 3 * n / 4; i++) {
        for (int j = 3 * m / 4 + 1; j < m; j++) {        
            check[i][j] += 3;
            v[i][j] += alpha * h[i][j];
        }
    }
    //4
    for (int i = 3 * n / 4+1; i < n ; i++) { 
        for (int j = 1; j < m; j++) {     
            check[i][j] += 4;
            v[i][j] += alpha * h[i][j];
        }
    }
    //5
    for (int i = n / 2; i <= 3 * n / 4; i++) {
        for (int j = 1; j < m / 4; j++) {
            check[i][j] += 5;
            v[i][j] += alpha * h[i][j];
        }
    }
    eps_max = 0;
    double beta = 0;
    S = 1;
    bool flag = false;
    while (!flag) {
        std::vector<std::vector<double> >  r(n + 1, std::vector<double>(m + 1));
      
        // 1
        for (int i = 1; i <= n / 4; i++) {
            for (int j = 1; j < 3 * m / 4; j++) {  
                check[i][j] += 1;
                r[i][j] = (a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
            }
        }
        //2
        for (int i = n / 4 + 1; i < n / 2; i++) {
            for (int j = 1; j < m; j++) { 
                check[i][j] += 2;
                r[i][j] = (a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
            }
        }
        //3
        for (int i = n / 2; i <= 3 * n / 4; i++) {
            for (int j = 3 * m / 4 + 1; j < m; j++) {     
                check[i][j] += 3;
                r[i][j] = (a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
            }
        }
        //4
        for (int i = 3 * n / 4+1; i < n ; i++) {
            for (int j = 1; j < m; j++) {   
                check[i][j] += 4;
                r[i][j] = (a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
            }
        }
        //5
        for (int i = n / 2; i <= 3 * n / 4; i++) {
            for (int j = 1; j < m / 4; j++) {
                check[i][j] += 5;
                r[i][j] = (a2 * v[i][j] + h2 * v[i - 1][j] + h2 * v[i + 1][j] + k2 * v[i][j - 1] + k2 * v[i][j + 1] - func(a + i * (b - a) / n, c + j * (d - c) / m));
            }
        }
     

        beta = calculateBetaCut(h2, k2, a2, r, h, func, n, m, (b - a) / n, (d - c) / m, a, c);

      /*
        // 1
        for (int i = 1; i <= n / 4; i++) {
            for (int j = 1; j < 3 * m / 4; j++) { 
                check[i][j] += 1;
                h[i][j] = beta * h[i][j] - r[i][j];
            }
        }
        //2
        for (int i = n / 4 + 1; i < n / 2; i++) {
            for (int j = 1; j < m; j++) {            
                check[i][j] += 2;
                h[i][j] = beta * h[i][j] - r[i][j];
            }
        }
        //3
        for (int i = n / 2; i <= 3 * n / 4; i++) {
            for (int j = 3 * m / 4 + 1; j < m; j++) {         
                check[i][j] += 3;
                h[i][j] = beta * h[i][j] - r[i][j];
            }
        }
        //4
        for (int i = 3 * n / 4+1; i < n ; i++) {
            for (int j = 1; j < m; j++) {       
                check[i][j] += 4;
                h[i][j] = beta * h[i][j] - r[i][j];
            }
        }
        //5
        for (int i = n / 2; i <= 3 * n / 4; i++) {
            for (int j = 1; j < m / 4; j++) {      
                check[i][j] += 5;
                h[i][j] = beta * h[i][j] - r[i][j];
            }
        }
        */
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m + 1; j++) {
                h[i][j] = beta * h[i][j] - r[i][j];
            }
        }
        
        alpha = calculateAlphaCut(h2, k2, a2, v, h, func, n, m, (b - a) / n, (d - c) / m, a, c);

        double save_component = 0;
        eps_max = 0;
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < m + 1; j++) {
                save_component = v[i][j];
                v[i][j] += alpha * h[i][j];
                save_component = abs(v[i][j] - save_component);
                if (save_component > eps_max)
                    eps_max = save_component;

            }
        }
     
        S += 1;
        if (S >= Nmax) {
            flag = true;
        }

    }
}

/*template <typename Function>

void solve(std::vector<std::vector<double>>& v, Function func, const int n, const int m, double a, double b, double c, double d, int Nmax, int& S, double& eps, double& eps_max, double& error_max) {
    double w = w_optimal(a, b, c, d, (b - a) / n, (d - c) / m);
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < m + 1; j++) {
            v[i].push_back(0);
        }
    }

    auto m1 = [](double x, double y) {
        return equalsZero(sin(M_PI * y) * sin(M_PI * y));
    };
    auto m2 = [](double x, double y) {
        return equalsZero(sin(M_PI * 2 * y) * sin(M_PI * 2 * y));
    };
    auto m3 = [](double x, double y) {
        return equalsZero(sin(M_PI * x) * sin(M_PI * x));
    };
    auto m4 = [](double x, double y) {
        return equalsZero(sin(M_PI * 2 * x) * sin(M_PI * 2 * x));
    };

    for (j = 0; j < m + 1; j++) {
        v[0][j] = m1(a, c + (d - c) / m * j);
    }
    for (j = 0; j < m + 1; j++) {
        v[n][j] = m2(b, c + (d - c) / m * j);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][0] = m3(a + i * (b - a) / n, c);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][m] = m4(a + i * (b - a) / n, d);
    }

    writeTable(n, m, v, a, b, c, d, 2);

    while (!flag) {
        eps_max = 0;
        for (j = 1; j < m; j++)
            for (i = 1; i < n; i++) {
                v_old = v[i][j];
                v_new = -w * (h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                double p = a + i * (b - a) / n;
                double e = c + j * (d - c) / m;
                double func_t = func(p, e);
                v_new = v_new + (1 - w) * a2 * v[i][j] + w * func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                if (eps_cur > eps_max) { eps_max = eps_cur; };
                v[i][j] = v_new;
            }
        if (S < 2) {
            writeTable(n, m, v, a, b, c, d, S * 2 + 4);
        }
        S = S + 1;
        if ((eps_max < eps) || (S >= Nmax)) { flag = true; }
    }
}*/




    

double presision(double x) {
    return x;
}

//int main() {
//    auto func = [](double x, double y) {
//        return abs(x * x - 2 * y);
//    };
//
//    auto u = [](double x, double y) {
//        return exp(sin(M_PI * x * y) * sin(M_PI * x * y));
//    };
//
//    int Nmax = 10000; // максимальное число итераций (не менее 1)
//    int S = 0; // счетчик итераций
//    int S_2 = 0; // счетчик итераций
//    double eps = 1e-8; // минимально допустимый прирост
//    double eps_max = 0; // текущее значение прироста
//    double eps_max_2 = 0; // текущее значение прироста
//    double eps_cur = 0; // для подсчета текущего значения прироста
//    double error_max = 0; // для подсчета текущего значения прироста
//    double accuracy = 1000; // для подсчета текущего значения прироста
//    double a2, k2, h2; // ненулевые элементы матрицы (-A)
//
//
//    const int n = 10, m = 10; //размерность сетки
//    writeHeader(n, m);
//
//    std::vector<std::vector<double>> v(n + 1); // искомый вектор v
//    std::vector<std::vector<double>> v_2(2 * n + 1); // искомый вектор v с половинным шагом
//    std::vector<double> r((n - 1) * (m - 1)); // невязка
//    double a = 0, b = 2, c = 0, d = 1; // границы области определения уравнения
//    double w = w_optimal(a, b, c, d, (b - a) / n, (d - c) / m);
//    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
//    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
//    a2 = -2 * (h2 + k2);
//    bool flag = false;
//
//    int i, j; //индексы
//    double v_old; // старое значение преобразуемой компоненты вектора v
//    double v_new; // новое значение преобразуемой компоненты вектора v
//
//    solve(v_2, func, 2 * n, 2 * m, a, b, c, d, Nmax, S, eps, eps_max, error_max);
//    S_2 = S;
//    S = 0;
//    eps_max_2 = eps_max;
//    eps_max = 0;
//    eps_cur = 0;
//    error_max = 0;
//    solve(v, func, n, m, a, b, c, d, Nmax, S, eps, eps_max, error_max);
//
//
//    for (j = 1; j < m; j++)
//        for (i = 1; i < n; i++) {
//            double k_1 = h2 * (v[i + 1][j] * ((i + 1 == n) ? 0 : 1)
//                + v[i - 1][j] * ((i - 1 == 0) ? 0 : 1));
//            double k_2 = k2 * (v[i][j + 1] * ((j + 1 == m) ? 0 : 1)
//                + v[i][j - 1] * ((j - 1 == 0) ? 0 : 1));
//            double k_3 = v[i][j] * a2;
//            v_new = (k_3 + k_1 + k_2);
//            double p = a + i * (b - a) / n;
//            double e = c + j * (d - c) / m;
//            double func_t = func(p, e);
//            double f_t = -(h2 * (v[i + 1][j] * ((i + 1 == n) ? 1 : 0)
//                + v[i - 1][j] * ((i - 1 == 0) ? 1 : 0)) +
//                k2 * (v[i][j + 1] * ((j + 1 == m) ? 1 : 0)
//                    + v[i][j - 1] * ((j - 1 == 0) ? 1 : 0))) + func_t;
//            r[(j - 1) * (m - 1) + (i - 1)] = v_new - f_t;
//        }
//
//    double r_norm = 0;
//    for (j = 0; j < r.size(); j++)
//        r_norm += r[j] * r[j];
//    r_norm = sqrt(r_norm);
//    std::ofstream outfile("test.dat");
//
//    for (j = 1; j < m; j++)
//        for (i = 1; i < n; i++) {
//            if (abs(v[i][j] - v_2[i][j]) < accuracy) {
//                accuracy = abs(v[i][j] - v_2[i][j]);
//            }
//        }
//
//    for (i = 0; i < n + 1; i++)
//        for (j = 0; j < m + 1; j++) {
//            v_new = v[i][j];
//            double p = a + i * (b - a) / n;
//            double e = c + j * (d - c) / m;
//            outfile << p << "\t" << e << "\t" << v_new << "\n";
//        }
//    //writeTable(n, m, vec_u, a, b, c, d, 8);
//    writeFinalTable(n, m, v, a, b, c, d, S, S_2, eps, eps_max, eps_max_2, error_max, r_norm, w, accuracy);
//    outfile.close();
//    system("python show_plot.py");
//}
auto funcTest = [](double x, double y) {
    return -((4 * M_PI * M_PI * y * y * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * cos(M_PI * x * y) * cos(M_PI * x * y) * sin(M_PI * x * y) * sin(M_PI * x * y) -
        2 * M_PI * M_PI * y * y * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * sin(M_PI * x * y) * sin(M_PI * x * y) +
        2 * M_PI * M_PI * y * y * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * cos(M_PI * x * y) * cos(M_PI * x * y)) +
        (4 * M_PI * M_PI * x * x * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * cos(M_PI * x * y) * cos(M_PI * x * y) * sin(M_PI * x * y) * sin(M_PI * x * y) -
            2 * M_PI * M_PI * x * x * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * sin(M_PI * x * y) * sin(M_PI * x * y) +
            2 * M_PI * M_PI * x * x * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * cos(M_PI * x * y) * cos(M_PI * x * y)));
};

auto funcDef = [](double x, double y) {
    return abs(x * x - 2 * y);
};

template <typename Function>
void solveMVR(std::vector<std::vector<double>>& v, Function func, const int n, const int m, double a, double b, double c, double d, int Nmax, int& S, double& eps, double& eps_max, double& error_max) {
    double w = w_optimal(a, b, c, d, (b - a) / n, (d - c) / m);
    int i, j; //èíäåêñû
    double a2, k2, h2; // íåíóëåâûå ýëåìåíòû ìàòðèöû (-A)
    double v_old; // ñòàðîå çíà÷åíèå ïðåîáðàçóåìîé êîìïîíåíòû âåêòîðà v
    double v_new; // íîâîå çíà÷åíèå ïðåîáðàçóåìîé êîìïîíåíòû âåêòîðà v
    bool flag = false; // óñëîâèå îñòàíîâêè
    double eps_cur = 0; // äëÿ ïîäñ÷åòà òåêóùåãî çíà÷åíèÿ ïðèðîñòà
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < m + 1; j++) {
            v[i].push_back(0);
        }
    }

    auto m1 = [](double x, double y) {
        return equalsZero(sin(M_PI * y) * sin(M_PI * y));
    };
    auto m2 = [](double x, double y) {
        return equalsZero(sin(M_PI * 2 * y) * sin(M_PI * 2 * y));
    };
    auto m3 = [](double x, double y) {
        return equalsZero(sin(M_PI * x) * sin(M_PI * x));
    };
    auto m4 = [](double x, double y) {
        return equalsZero(sin(M_PI * 2 * x) * sin(M_PI * 2 * x));
    };

    for (j = 0; j < m + 1; j++) {
        v[0][j] = m1(a, c + (d - c) / m * j);
    }
    for (j = 0; j < m + 1; j++) {
        v[n][j] = m2(b, c + (d - c) / m * j);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][0] = m3(a + i * (b - a) / n, c);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][m] = m4(a + i * (b - a) / n, d);
    }

    writeTable(n, m, v, a, b, c, d, 2);

    while (!flag) {
        eps_max = 0;
        for (j = 1; j < m; j++)
            for (i = 1; i < n; i++) {
                v_old = v[i][j];
                v_new = -w * (h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                double p = a + i * (b - a) / n;
                double e = c + j * (d - c) / m;
                double func_t = func(p, e);
                v_new = v_new + (1 - w) * a2 * v[i][j] + w * func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                if (eps_cur > eps_max) { eps_max = eps_cur; };
                v[i][j] = v_new;
            }
        if (S < 2) {
            writeTable(n, m, v, a, b, c, d, S * 2 + 4);
        }
        S = S + 1;
        if ((eps_max < eps) || (S >= Nmax)) { flag = true; }
    }
}


struct func {
    double operator()(double x, double y) {
        return abs(x * x - 2 * y);
    }

    void operator()() {
        system("python show_plot.py");
    }
};
