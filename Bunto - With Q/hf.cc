#include"hf.h"


double integrate(std::function<double(double)> f, double a, double b, double acc, double eps, double f2, double f3) // NaN indicates first call
{
    double h = b - a;
    if(std::isnan(f2)){f2 = f(a+2*h/6); f3 = f(a+4*h/6);} // first call, no points to reuse
    double f1 = f(a+h/6.), f4 = f(a+5.*h/6.);
    double Q = (2.*f1+f2+f3+2.*f4)/6.*(b - a); // higher order rule
    double q = (  f1+f2+f3+  f4)/4.*(b - a); // lower order rule
    double err = std::abs(Q - q);
    if (err <= acc + eps*std::abs(Q)){return Q;}
    else {return integrate(f, a, (a + b)/2., acc/std::sqrt(2), eps, f1, f2) + integrate(f, (a + b)/2., b, acc/std::sqrt(2), eps, f3, f4);};
};


double integrate_CC(std::function<double(double)> f, double a, double b, double acc, double eps){
    if (std::isinf(a) && a < 0 && std::isinf(b) && b > 0){ // eq. 62
        auto F = [=](double t){
            if (t == 0){return 0.;}
            double value = (1. - t)/t;
            return (f(value) + f(-value))/t/t;
        };
        return integrate_CC(F, 0., 1., acc, eps);
    }
    if (std::isinf(a) && a < 0){ // eq. 66
        auto F = [=](double t){
            if (t == 0){return 0.;}
            double value = (1. - t)/t;
            return f(b - value)/t/t;
        };
        return integrate_CC(F, 0., 1., acc, eps);
    }
    if (std::isinf(b) && b > 0){ // eq. 64
        auto F = [=](double t){
            if (t == 0){return 0.;}
            double value = (1. - t)/t;
            return f(a + value)/t/t;
        };
        return integrate_CC(F, 0., 1., acc, eps);
    }
    auto F = [=](double theta){
        double value = (a + b)/2. + (b - a)/2.*std::cos(theta);
        return f(value)*std::sin(theta)*(b - a)/2.;
    };
    return integrate(F, 0., pi, acc, eps);
};


std::complex<double> integrate_CC_C(const std::function<std::complex<double>(double)>& f, double a, double b, int N){
    if (std::isinf(a) && a < 0 && std::isinf(b) && b > 0){ // eq. 62
        auto F = [=](double t){
            if (t == 0){return std::complex<double>{};}
            double value = (1. - t)/t;
            return (f(value) + f(-value))/t/t;
        };
        return integrate_CC_C(F, 0., 1., N);
    }
    if (std::isinf(a) && a < 0){ // eq. 66
        auto F = [=](double t){
            if (t == 0){return std::complex<double>{};}
            double value = (1. - t)/t;
            return f(b - value)/t/t;
        };
        return integrate_CC_C(F, 0., 1., N);
    }
    if (std::isinf(b) && b > 0){ // eq. 64
        auto F = [=](double t){
            if (t == 0){return std::complex<double>{};}
            double value = (1. - t)/t;
            return f(a + value)/t/t;
        };
        return integrate_CC_C(F, 0., 1., N);
    }
    auto F = [=](double theta){
        double value = (a + b)/2. + (b - a)/2.*std::cos(theta);
        return f(value)*std::sin(theta)*(b - a)/2.;
    };
    return integrate_trapezoid(F, 0., pi, N);
};


std::complex<double> integrate_trapezoid(std::function<std::complex<double>(double)> f, double a, double b, int N) // NaN indicates first call
{
    double h =(b - a)/(N - 1);
    std::vector<double> xs = linspace(a, b, N);
    std::complex<double> f1, f2, sum = f(xs[0])/2. + f(xs[N - 1])/2.;
    for (int i = 1; i < N - 1; ++i){
        sum += f(xs[i]);
    }
    return sum*h;
};


std::complex<double> C_integrate_CC(const std::function<std::complex<double>(std::complex<double>)>& f, double a1, double b1, double a2, double b2, int N, int M){
    auto inner = [&] (double x){
        auto F = [&] (double y){return f(std::complex<double>{x, y});};
        return integrate_CC_C(F, a2, b2, M);
    };
    return integrate_CC_C(inner, a1, b1, N);
}