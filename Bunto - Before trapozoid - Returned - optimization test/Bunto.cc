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
    double cosh_r = std::cosh(r), norm_alpha_0 = std::norm(alpha_0), tanh_r = std::tanh(r), prefactor = 1/(pi*cosh_r);
    double su = 1/std::sqrt(1 + tanh_r), sv = 1/std::sqrt(1 - tanh_r);
    alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
    std::complex<double> phase = std::polar(1., std::arg(z)), conj_phase = std::conj(phase), conj_alpha_0 = std::conj(alpha_0);
    double x0 = std::abs(alpha_0) * std::exp(-r), L = 5, M = 3.25; //L = 2.8, M = 2.6; // L = 5, M = 3.25;

    std::function<double(std::complex<double>)> Q = [&](std::complex<double> beta){
        std::complex<double> conj_beta = std::conj(beta);
        return std::real(prefactor*std::exp(-(norm_alpha_0 + std::norm(beta)) + (conj_beta*alpha_0 + beta*conj_alpha_0)/cosh_r - ((phase*(conj_beta*conj_beta - conj_alpha_0*conj_alpha_0) + conj_phase*(beta*beta - alpha_0*alpha_0))*tanh_r/2.)));
    };
    std::function<std::complex<double>(std::complex<double>)> f_a0_squared = [&](std::complex<double> z){
       
        std::complex<double> beta = x0 + su*std::real(z) + std::complex<double>{0., 1.}* sv*std::imag(z);
        return Q(beta)*beta*beta*su*sv;
    };
    std::function<std::complex<double>(std::complex<double>)> f_a0_dagger_a0 = [&](std::complex<double> z){
        std::complex<double> beta = x0 + su*std::real(z) + std::complex<double>{0., 1.}* sv*std::imag(z);
        return Q(beta)*(std::norm(beta) - 1)*su*sv;
    };

    double acc = 1e-2, acc2 = acc;
    
    a0_squared = C_integrate_CC(f_a0_squared, x0 - M, x0 + M, -L, L, acc, acc);
    // test = C_integrate_CC(f_a0_dagger_a0, x0 - M, x0 + M, -L, L, acc2);
    // a0_dagger_a0 = std::real(test);
    a0_dagger_a0 = std::norm(alpha_0) + std::sinh(r)*std::sinh(r);
}