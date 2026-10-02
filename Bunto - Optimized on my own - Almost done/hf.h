#pragma once
#include<functional>
#include<limits>
#include<tuple>
#include<complex>
#include<cmath>
#include<vector>
#include<iostream>

template <typename T> std::vector<T> linspace(T nmin, T nmax, int N){
    std::vector<T> v(N);
    T grad = (nmax - nmin)/(N - 1);
    for (int i = 0; i < N; ++i){
        v[i] = i*grad + nmin;
    }
    return v;
}
template <typename T> std::vector<T> arange(T nmin, T nmax, T N){
    int M = (int)((nmax - nmin)/N);
    if (M <= 0){throw std::runtime_error("Cannot make arange with no elements");}
    std::vector<T> v(M);
    for (int i = 0; i < M; ++i){
        v[i] = nmin + i*N;
    }
    return v;
}

inline double NaN = std::numeric_limits<double>::quiet_NaN();
inline const double pi = 3.1415926535897932384626433;
inline const double inv_sqrt2 = 0.7071067811865475;
double integrate(std::function<double(double)> f, double a, double b, double acc = 0.001, double eps = 0.001, double f2 = NaN, double f3 = NaN); // NaN indicates first call
double integrate_CC(std::function<double(double)> f, double a, double b, double acc = 0.001, double eps = 0.001);
// std::complex<double> integrate_C(std::function<std::complex<double>(double)> f, double a, double b, double acc = 0.001, double eps = 0.001, std::complex<double> f2 = std::complex{0., 0.}, std::complex<double> f3 = std::complex{0., 0.}, bool is_first_call = true); // NaN indicates first call
std::complex<double> integrate_CC_C(const std::function<std::complex<double>(double)>& f, double a, double b, int N = 5e4);





template <typename T> std::complex<double> integrate_C(const T& f, double a, double b, double acc = 1e-3, double eps = 1e-3, std::complex<double> f2 = {}, std::complex<double> f3 = {}, int N = 0) // NaN indicates first call
{
    double h = b - a;
    if(N == 0){f2 = f(a+2.*h/6.); f3 = f(a+4.*h/6.);} // first call, no points to reuse
    std::complex<double> f1 = f(a+h/6.), f4 = f(a+5.*h/6.);
    std::complex<double> Q = (2.*f1+f2+f3+2.*f4)/6.*(b - a); // higher order rule
    std::complex<double> q = (  f1+f2+f3+  f4)/4.*(b - a); // lower order rule
    double err = std::abs(Q - q);
    if (err <= acc + eps*std::abs(Q) or N == 25){return Q;}
    else {return integrate_C(f, a, (a + b)/2., acc*inv_sqrt2, eps, f1, f2, N + 1) + integrate_C(f, (a + b)/2., b, acc*inv_sqrt2, eps, f3, f4, N + 1);};
};



std::complex<double> integrate_trapezoid(std::function<std::complex<double>(double)> f, double a, double b, int N = 5e4); // NaN indicates first call







std::complex<double> C_integrate_CC(const std::function<std::complex<double>(std::complex<double>)>& f, double a1, double b1, double a2, double b2, int N = 5e4, int M = 5e4);
