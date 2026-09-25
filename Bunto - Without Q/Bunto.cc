#include"Bunto.h"
#include<cmath>
#include<iostream>
#include<complex>


// double sim::A(double E){
//     double norm_Me = std::norm(M(E, 1)), norm_Ma = std::norm(M(E, 0));
//     return 1/4.*f0*f0*a0_dagger_a0 * (norm_Me*Fe*Fe + norm_Ma*Fa*Fa);
// } // Amplitude
// double sim::C(double E){
//     double norm = std::abs(xi(E));
//     return 2*norm*mu/(norm*norm*mu*mu + 1)*std::abs(a0_squared)/a0_dagger_a0;
// } // Contrast

// double sim::theta(double E){return -std::arg(xi(E)) + std::arg(a0_squared);} // Phase
// std::complex<double> sim::xi(double E){return fe*M(E, 1)/(fa*M(E, 0));}
// std::complex<double> sim::M(double E, int beta){
//     if(E == 1.892 && beta == 0){return std::polar(86.8374, 0.3692);}
//     if(E == 1.892 && beta == 1){return std::polar(56.7437, -0.3165);}
//     if(E == 4.990 && beta == 0){return std::polar(79.2683, 0.2293);}
//     if(E == 4.990 && beta == 1){return std::polar(55.1106, -0.2117);}
//     if(E == 8.088 && beta == 0){return std::polar(67.3156, 0.1585);}
//     if(E == 8.088 && beta == 1){return std::polar(48.4818, -0.1542);}
//     if(E == 11.19 && beta == 0){return std::polar(55.9270, 0.1192);}
//     if(E == 11.19 && beta == 1){return std::polar(41.6357, -0.1188);}
//     if(E == 14.29 && beta == 0){return std::polar(46.4232, 0.0944);}
//     if(E == 14.29 && beta == 1){return std::polar(35.4227, -0.0953);}
//     if(E == 17.38 && beta == 0){return std::polar(38.9820, 0.0769);}
//     if(E == 17.38 && beta == 1){return std::polar(30.1262, -0.0788);}
//     if(E == 20.48 && beta == 0){return std::polar(32.7244, 0.0647);}
//     if(E == 20.48 && beta == 1){return std::polar(25.7704, -0.0666);}
//     if(E == 23.58 && beta == 0){return std::polar(27.5185, 0.0560);}
//     if(E == 23.58 && beta == 1){return std::polar(22.1855, -0.0571);}
//     throw std::runtime_error("Not an allowed energy");
// } // Second order matrix element
// double sim::f(double omega){return std::sqrt(omega/(2*epsilon_0*V));} // Vacuum electric field amplitude
// double sim::I(double E, double tau){return A(E)*(1 + C(E)*std::cos(2*omega_0*tau + theta(E)));}

// void sim::update(std::complex<double> z){
//     // zeta = z; 
//     // r = std::abs(z);
//     // alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
//     // a0_squared = alpha_0*alpha_0 - std::polar(1., std::arg(z))/2.*std::sinh(2*r);
//     // a0_dagger_a0 = std::norm(alpha_0) + std::sinh(r)*std::sinh(r);


//     zeta = z; 
//     r = std::abs(z);
//     alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
//     std::complex<double> phase = std::polar(1., std::arg(z));
//     std::function<std::complex<double>(std::complex<double>)> Q = [&](std::complex<double> beta){
//         std::complex<double> conj_beta = std::conj(beta), conj_alpha_0 = std::conj(alpha_0);
//         return std::exp(-(std::norm(alpha_0) + std::norm(beta)) - (conj_beta*alpha_0 + beta*conj_alpha_0)/std::cosh(r) - (phase*(conj_beta*conj_beta - conj_alpha_0*conj_alpha_0) + std::conj(phase)*(beta*beta - alpha_0*alpha_0)*std::tanh(r)/2.));
//     };
//     std::function<std::complex<double>(std::complex<double>)> f_a0_squared = [&](std::complex<double> x){return x*x;};
//     std::function<std::complex<double>(std::complex<double>)> f_a0_dagger_a0 = [&](std::complex<double> x){return std::norm(x);};
//     // std::complex<double> f_a0_dagger_a0 = [&](std::complex<double>){return std::norm(alpha_0);};

    
//     // double inf = std::numeric_limits<double>::infinity();
//     a0_squared = C_integrate_CC(f_a0_squared, -inf, inf);
//     // a0_dagger_a0 = C_integrate_CC(f_a0_dagger_a0, -inf, inf);
//     a0_dagger_a0 = std::norm(alpha_0) + std::sinh(r)*std::sinh(r);
// }



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
    alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
    a0_squared = alpha_0*alpha_0 - std::polar(1., std::arg(z))/2.*std::sinh(2*r);
    a0_dagger_a0 = std::norm(alpha_0) + std::sinh(r)*std::sinh(r);
}












// #include"Bunto.h"
// #include<cmath>
// #include<iostream>
// #include<complex>

// double pi = 3.14159265358979323846;

// // std::sqrt(omega/(2*epsilon_0*V))
// // Note that norm(z) returns norm squared of z, while abs(z) returns the norm.
// double n = 1e16;
// double omega_0 = 0.057, epsilon_0 = 1, V = omega_0/(2*epsilon_0*1e-24), varepsilon_0 = 0.5, r = std::asinh(std::sqrt(varepsilon_0*n)); // roughly 800 nm
// double omega_a = 7*omega_0, omega_e = omega_a + 2*omega_0, fa = f(omega_a), fe = f(omega_e), f0 = f(omega_0);


// // double Fe = 2*fe*std::sqrt(n), Fa = 2*fa*std::sqrt(n);
// double Fe = 1.7e-5, Fa = 1.7e-5;
// double mu = Fe*fa/(Fa*fe);

// std::vector<double> results(double E, std::complex(zeta)){

// }
// std::complex<double> zeta = std::polar(r, pi), alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0); // Modulus and arg
// std::complex<double> a0_squared = alpha_0*alpha_0 - zeta/r/2.*std::sinh(2*r);
// double a0_dagger_a0 = std::norm(alpha_0) + std::sinh(r)*std::sinh(r);
// double A(double E){
//     double norm_Me = std::norm(M(E, 1)), norm_Ma = std::norm(M(E, 0));
//     return 1/4.*f0*f0*a0_dagger_a0 * (norm_Me*Fe*Fe + norm_Ma*Fa*Fa)*0.7729;
// } // Amplitude
// double C(double E){
//     double norm = std::abs(xi(E));
//     return 2*norm*mu/(norm*norm*mu*mu + 1)*std::abs(a0_squared)/a0_dagger_a0;
// } // Contrast
// double theta(double E){return -std::arg(xi(E)) + std::arg(a0_squared);} // Phase
// std::complex<double> xi(double E){return fe*M(E, 1)/(fa*M(E, 0));}
// std::complex<double> M(double E, int beta){
//     if(E == 1.892 && beta == 0){return std::polar(86.8374, 0.3692);}
//     if(E == 1.892 && beta == 1){return std::polar(56.7437, -0.3165);}
//     if(E == 4.990 && beta == 0){return std::polar(79.2683, 0.2293);}
//     if(E == 4.990 && beta == 1){return std::polar(55.1106, -0.2117);}
//     if(E == 8.088 && beta == 0){return std::polar(67.3156, 0.1585);}
//     if(E == 8.088 && beta == 1){return std::polar(48.4818, -0.1542);}
//     if(E == 11.19 && beta == 0){return std::polar(55.9270, 0.1192);}
//     if(E == 11.19 && beta == 1){return std::polar(41.6357, -0.1188);}
//     if(E == 14.29 && beta == 0){return std::polar(46.4232, 0.0944);}
//     if(E == 14.29 && beta == 1){return std::polar(35.4227, -0.0953);}
//     if(E == 17.38 && beta == 0){return std::polar(38.9820, 0.769);}
//     if(E == 17.38 && beta == 1){return std::polar(30.1262, -0.0788);}
//     if(E == 20.48 && beta == 0){return std::polar(32.7244, 0.0647);}
//     if(E == 20.48 && beta == 1){return std::polar(25.7704, -0.0666);}
//     if(E == 23.58 && beta == 0){return std::polar(27.5185, 0.0560);}
//     if(E == 23.58 && beta == 1){return std::polar(22.1855, -0.0571);}
//     throw std::runtime_error("Not an allowed energy");
// } // Second order matrix element
// double f(double omega){return std::sqrt(omega/(2*epsilon_0*V));} // Vacuum electric field amplitude
// double I(double E, double tau){return A(E)*(1 + C(E)*std::cos(2*omega_0*tau + theta(E)));}