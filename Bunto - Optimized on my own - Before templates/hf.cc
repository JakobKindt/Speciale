#include"hf.h"

inline double pi = 3.1415926535897932384626433;
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






std::complex<double> integrate_C(std::function<std::complex<double>(double)> f, double a, double b, double acc, double eps, std::complex<double> f2, std::complex<double> f3, bool is_first_call) // NaN indicates first call
{
    double h = b - a;
    if(is_first_call){f2 = f(a+2*h/6.); f3 = f(a+4.*h/6.);} // first call, no points to reuse
    std::complex<double> f1 = f(a+h/6.), f4 = f(a+5.*h/6.);
    std::complex<double> Q = (2.*f1+f2+f3+2.*f4)/6.*(b - a); // higher order rule
    std::complex<double> q = (  f1+f2+f3+  f4)/4.*(b - a); // lower order rule
    double err = std::abs(Q - q);
    if (err <= acc + eps*std::abs(Q)){return Q;}
    else {return integrate_C(f, a, (a + b)/2., acc/std::sqrt(2), eps, f1, f2, false) + integrate_C(f, (a + b)/2., b, acc/std::sqrt(2), eps, f3, f4, false);};
};


std::complex<double> integrate_CC_C(std::function<std::complex<double>(double)> f, double a, double b, double acc, double eps){
    if (std::isinf(a) && a < 0 && std::isinf(b) && b > 0){ // eq. 62
        auto F = [=](double t){
            if (t == 0){return std::complex{0., 0.};}
            double value = (1. - t)/t;
            return (f(value) + f(-value))/t/t;
        };
        return integrate_CC_C(F, 0., 1., acc, eps);
    }
    if (std::isinf(a) && a < 0){ // eq. 66
        auto F = [=](double t){
            if (t == 0){return std::complex{0., 0.};}
            double value = (1. - t)/t;
            return f(b - value)/t/t;
        };
        return integrate_CC_C(F, 0., 1., acc, eps);
    }
    if (std::isinf(b) && b > 0){ // eq. 64
        auto F = [=](double t){
            if (t == 0){return std::complex{0., 0.};}
            double value = (1. - t)/t;
            return f(a + value)/t/t;
        };
        return integrate_CC_C(F, 0., 1., acc, eps);
    }
    auto F = [=](double theta){
        double value = (a + b)/2. + (b - a)/2.*std::cos(theta);
        return f(value)*std::sin(theta)*(b - a)/2.;
    };
    return integrate_C(F, 0., pi, acc, eps);
};












std::complex<double> C_integrate_CC(std::function<std::complex<double>(std::complex<double>)> f, double a1, double b1, double a2, double b2, double acc, double eps){
    auto inner = [&] (double x){
        auto F = [&] (double y){return f(std::complex{x, y});};
        return integrate_CC_C(F, a2, b2, acc/std::sqrt(2), eps);
    };
    return integrate_CC_C(inner, a1, b1, acc/std::sqrt(2), eps);
}




// std::complex<double> C_integrate_CC(std::function<std::complex<double>(std::complex<double>)> f, double a1, double b1, double a2, double b2, double acc, double eps){
//     double re, im;
//     auto inner_re_f = [&] (double x){
//         auto F = [&] (double y){return std::real(f(std::complex{x, y}));};
//         return integrate_CC(F, a2, b2, acc/std::sqrt(2), eps);
//     };

//     auto inner_im_f = [&] (double x){
//         auto F = [&] (double y){return std::imag(f(std::complex{x, y}));};
//         return integrate_CC(F, a2, b2, acc/std::sqrt(2), eps);
//     };
//     re = integrate_CC(inner_re_f, a1, b1, acc/std::sqrt(2), eps);
//     im = integrate_CC(inner_im_f, a1, b1, acc/std::sqrt(2), eps);
//     return std::complex{re, im};
// }