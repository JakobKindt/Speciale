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
double integrate(std::function<double(double)> f, double a, double b, int N = 1e4); // NaN indicates first call
// double integrate(std::function<double(double)> f, double a, double b, double acc = 0.001, double eps = 0.001); // NaN indicates first call
double integrate_CC(std::function<double(double)> f, double a, double b, int N = 1e4);
// double integrate_CC(std::function<double(std::complex<double>)> f, double a, double b, double acc = 0.001, double eps = 0.001);
std::complex<double> C_integrate_CC(std::function<std::complex<double>(std::complex<double>)> f, double a1, double b1, double a2, double b2, int N = 1e4, int M = 1e4);
