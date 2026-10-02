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
inline double pi = 3.1415926535897932384626433;
template <typename T, typename R> R integrate(T& f, double a, double b, double acc = 0.001, double eps = 0.001, R f2 = R{}, R f3 = R{}, bool is_first_call = true) // NaN indicates first call
{
    double h = b - a;

    if(is_first_call){f2 = f(a+2*h/6); f3 = f(a+4*h/6);} // first call, no points to reuse
    R f1 = f(a+h/6), f4 = f(a+5*h/6);
    R Q = (2.*f1+f2+f3+2.*f4)/6.*(b - a); // higher order rule
    R q = (  f1+f2+f3+  f4)/4.*(b - a); // lower order rule
    double err = std::abs(Q - q);
    if (err <= acc + eps*std::abs(Q)){return Q;}
    else {return integrate<T, R>(f, a, (a + b)/2., acc/std::sqrt(2), eps, f1, f2, false) + integrate<T, R>(f, (a + b)/2., b, acc/std::sqrt(2), eps, f3, f4, false);};
};


template <typename T, typename R> R integrate_CC(T& f, double a, double b, double acc = 0.001, double eps = 0.001){
    if (std::isinf(a) && a < 0 && std::isinf(b) && b > 0){ // eq. 62
        auto F = [=](double t){
            if (t == 0){return R{};} // R{} = 0, even for complex case.
            double value = (1. - t)/t;
            return (f(value) + f(-value))/t/t;
        };
        return integrate_CC<decltype(F), R>(F, 0., 1., acc, eps);
    }
    if (std::isinf(a) && a < 0){ // eq. 66
        auto F = [=](double t){
            if (t == 0){return R{};}
            double value = (1. - t)/t;
            return f(b - value)/t/t;
        };
        return integrate_CC<decltype(F), R>(F, 0., 1., acc, eps);
    }
    if (std::isinf(b) && b > 0){ // eq. 64
        auto F = [=](double t){
            if (t == 0){return R{};}
            double value = (1 - t)/t;
            return f(a + value)/t/t;
        };
        return integrate_CC<decltype(F), R>(F, 0., 1., acc, eps);
    }
    auto F = [=](double theta){
        double value = (a + b)/2 + (b - a)/2*std::cos(theta);
        return f(value)*std::sin(theta)*(b - a)/2.;
    };
    return integrate<decltype(F), R>(F, 0., pi, acc, eps);
};


std::complex<double> C_integrate_CC(std::function<std::complex<double>(std::complex<double>)>& f, double a1, double b1, double a2, double b2, double acc = 0.001, double eps = 0.001){
    auto inner = [&] (double x){
        auto F = [&] (double y){return f(std::complex{x, y});};
        return integrate_CC<decltype(F), std::complex<double>>(F, a2, b2, acc/std::sqrt(2), eps);
    };
    return integrate_CC<decltype(inner), std::complex<double>>(inner, a1, b1, acc/std::sqrt(2), eps);
}





















// double integrate(std::function<double(double)> f, double a, double b, double acc = 0.001, double eps = 0.001, double f2 = NaN, double f3 = NaN); // NaN indicates first call
// double integrate_CC(std::function<double(double)> f, double a, double b, double acc = 0.001, double eps = 0.001);
// std::tuple<double, double> integrate_with_err(std::function<double(double)> f, double a, double b, double acc = 0.001, double eps = 0.001, double f2 = NaN, double f3 = NaN);
// // double integrate_CC(std::function<double(std::complex<double>)> f, double a, double b, double acc = 0.001, double eps = 0.001);
// std::complex<double> C_integrate_CC(std::function<std::complex<double>(std::complex<double>)> f, double a1, double b1, double a2, double b2, double acc = 0.001, double eps = 0.001);
