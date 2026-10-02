#include"Bunto.h"
#include<cmath>
#include<iostream>
#include<complex>


double sim::A(double E){
    double norm_Me = std::norm(M(E, 1)), norm_Ma = std::norm(M(E, 0));
    return 1/4.*f0*f0*a0_dagger_a0 * (norm_Me*Fe*Fe + norm_Ma*Fa*Fa);
} // Amplitude
double sim::C(double E){
    double norm = std::abs(xi(E));
    return 2*norm*mu/(norm*norm*mu*mu + 1)*std::abs(a0_squared)/a0_dagger_a0;
} // Contrast

double sim::theta(double E){return -std::arg(xi(E)) + std::arg(a0_squared);} // Phase
std::complex<double> sim::xi(double E){return fe*M(E, 1)/(fa*M(E, 0));}
std::complex<double> sim::M(double E, int beta){
    if(E == 1.892 && beta == 0){return std::polar(86.8374, 0.3692);}
    if(E == 1.892 && beta == 1){return std::polar(56.7437, -0.3165);}
    if(E == 4.990 && beta == 0){return std::polar(79.2683, 0.2293);}
    if(E == 4.990 && beta == 1){return std::polar(55.1106, -0.2117);}
    if(E == 8.088 && beta == 0){return std::polar(67.3156, 0.1585);}
    if(E == 8.088 && beta == 1){return std::polar(48.4818, -0.1542);}
    if(E == 11.19 && beta == 0){return std::polar(55.9270, 0.1192);}
    if(E == 11.19 && beta == 1){return std::polar(41.6357, -0.1188);}
    if(E == 14.29 && beta == 0){return std::polar(46.4232, 0.0944);}
    if(E == 14.29 && beta == 1){return std::polar(35.4227, -0.0953);}
    if(E == 17.38 && beta == 0){return std::polar(38.9820, 0.0769);}
    if(E == 17.38 && beta == 1){return std::polar(30.1262, -0.0788);}
    if(E == 20.48 && beta == 0){return std::polar(32.7244, 0.0647);}
    if(E == 20.48 && beta == 1){return std::polar(25.7704, -0.0666);}
    if(E == 23.58 && beta == 0){return std::polar(27.5185, 0.0560);}
    if(E == 23.58 && beta == 1){return std::polar(22.1855, -0.0571);}
    throw std::runtime_error("Not an allowed energy");
} // Second order matrix element
double sim::f(double omega){return std::sqrt(omega/(2*epsilon_0*V));} // Vacuum electric field amplitude
double sim::I(double E, double tau){return A(E)*(1 + C(E)*std::cos(2*omega_0*tau + theta(E)));}



void sim::update(std::complex<double> z){
    zeta = z; 
    r = std::abs(z); 
    varepsilon_0 = 0.5;
    if (r == 0){varepsilon_0 = 0.;}
    omega_0 = 0.057; epsilon_0 = 1; V = omega_0/(2*epsilon_0*1e-24); omega_a = 7*omega_0; omega_e = omega_a + 2*omega_0; fa = f(omega_a); fe = f(omega_e); f0 = f(omega_0);
    double theta = std::arg(z);
    alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
    
    // double su = 1/std::sqrt(1 + tanh_r), sv = 1/std::sqrt(1 - tanh_r); // susv = su*sv = cosh(r), so prefactor*cosh(r) = 1/pi which is more nummerical stable;
    double su = std::sqrt(0.5*(std::exp(-2.0*r) + 1.0)), sv = std::sqrt(0.5*(std::exp(2.0*r) + 1.0)); // ChatGPT. Verify later
    std::complex<double> re, im;
    double cos_2 = std::cos(theta/2), sin_2 = std::sin(theta/2);
    
    double L = 10, M = 10;
    auto Q = [&](std::complex<double> beta){
        return std::exp(-std::norm(beta))/pi;
    };
    auto f_a0_squared = [&](std::complex<double> z){
        re = std::real(z); im = std::imag(z);
        std::complex<double> beta = alpha_0 + (su*cos_2*re - sv*sin_2*im) + std::complex<double>{0., 1.}*(su*sin_2*re + sv*cos_2*im);
        return Q(z)*beta*beta;
    };
    auto f_a0_dagger_a0 = [&](std::complex<double> z){
        re = std::real(z); im = std::imag(z);
        std::complex<double> beta = alpha_0 + (su*cos_2*re - sv*sin_2*im) + std::complex<double>{0., 1.}*(su*sin_2*re + sv*cos_2*im);
        return Q(z)*(std::norm(beta) - 1);
    };

    double acc = 1e3, acc2 = acc*1e1;
    
    a0_dagger_a0 = std::real(C_integrate_CC(f_a0_dagger_a0, -M, M, -L, L, acc, acc));
    a0_squared = C_integrate_CC(f_a0_squared, -M, M, -L, L, acc2, acc2);
}