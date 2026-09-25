#include<complex>
#include"hf.h"
#include<functional>

inline double inf = std::numeric_limits<double>::infinity();


class sim{
public:
    double A(double E); // Amplitude
    double C(double E); // Contrast
    double theta(double E); // Phase
    std::complex<double> xi(double E);
    double f(double omega); // Vacuum electric field amplitude
    double I(double E, double tau = 0);
    double acc = 1e-10;
    std::complex<double> zeta, test;
    sim(std::complex<double> z){
        zeta = z; 
        r = std::abs(z);
        alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
        std::complex<double> phase = std::polar(1., std::arg(z));
        // std::function<std::complex<double>(std::complex<double>)> Q = [&](std::complex<double> beta){
        //     std::complex<double> conj_beta = std::conj(beta), conj_alpha_0 = std::conj(alpha_0);
        //     return 1/(pi*std::cosh(r))*std::exp(-(std::norm(alpha_0) + std::norm(beta)) + (conj_beta*alpha_0 + beta*conj_alpha_0)/std::cosh(r) - ((phase*(conj_beta*conj_beta - conj_alpha_0*conj_alpha_0) + std::conj(phase)*(beta*beta - alpha_0*alpha_0))*std::tanh(r)/2.));
        // };
        std::function<double(std::complex<double>)> Q = [&](std::complex<double> beta){
            std::complex<double> conj_beta = std::conj(beta), conj_alpha_0 = std::conj(alpha_0);
            return std::real(1/(pi*std::cosh(r))*std::exp(-(std::norm(alpha_0) + std::norm(beta)) + (conj_beta*alpha_0 + beta*conj_alpha_0)/std::cosh(r) - ((phase*(conj_beta*conj_beta - conj_alpha_0*conj_alpha_0) + std::conj(phase)*(beta*beta - alpha_0*alpha_0))*std::tanh(r)/2.)));
        };
        // std::complex<double> f_a0_squared = [&](std::complex<double>){return alpha_0*alpha_0;};
        // std::complex<double> f_a0_dagger_a0 = [&](std::complex<double>){return std::norm(alpha_0);};
        std::function<std::complex<double>(std::complex<double>)> f_a0_squared = [&](std::complex<double> x){return Q(x)*x*x;};
        std::function<std::complex<double>(std::complex<double>)> f_a0_dagger_a0 = [&](std::complex<double> x){return Q(x)*std::complex{std::norm(x) - 1., 0.};};

        // double inf = std::numeric_limits<double>::infinity();
        // double inf = std::numeric_limits<double>::infinity();
        a0_squared = C_integrate_CC(f_a0_squared, -inf, inf, acc, acc);
        test = C_integrate_CC(f_a0_dagger_a0, -inf, inf, acc, acc);
        a0_dagger_a0 = std::real(test);
        
        // a0_dagger_a0 = std::norm(alpha_0) + std::sinh(r)*std::sinh(r);
    }
    void update(std::complex<double> zeta);

// private:
    std::function<double(std::complex<double>)> Q;
    double pi = 3.14159265358979323846;
    double n = 1e16;
    double omega_0 = 0.057, epsilon_0 = 1, V = omega_0/(2*epsilon_0*1e-24), varepsilon_0 = 0.5; // roughly 800 nm
    double omega_a = 7*omega_0, omega_e = omega_a + 2*omega_0, fa = f(omega_a), fe = f(omega_e), f0 = f(omega_0);
    double Fe = 1.7e-5, Fa = 1.7e-5;
    double mu = Fe*fa/(Fa*fe);
    std::complex<double> alpha_0, a0_squared;
    double r, a0_dagger_a0;
    
    
    std::complex<double> M(double E, int beta); // Second order matrix element
};









// class sim{
// public:
//     double A(double E); // Amplitude
//     double C(double E); // Contrast
//     double theta(double E); // Phase
//     std::complex<double> xi(double E);
//     double f(double omega); // Vacuum electric field amplitude
//     double I(double E, double tau = 0);
//     std::complex<double> zeta;
//     sim(std::complex<double> z){
//         zeta = z; 
//         r = std::abs(z);
//         alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
//         a0_squared = alpha_0*alpha_0 - std::polar(1., std::arg(z))/2.*std::sinh(2*r);
//         a0_dagger_a0 = std::norm(alpha_0) + std::sinh(r)*std::sinh(r);
//     }
//     void update(std::complex<double> zeta);

// // private:

//     double pi = 3.14159265358979323846;
//     double n = 1e16;
//     double omega_0 = 0.057, epsilon_0 = 1, V = omega_0/(2*epsilon_0*1e-24), varepsilon_0 = 0.5; // roughly 800 nm
//     double omega_a = 7*omega_0, omega_e = omega_a + 2*omega_0, fa = f(omega_a), fe = f(omega_e), f0 = f(omega_0);
//     double Fe = 1.7e-5, Fa = 1.7e-5;
//     double mu = Fe*fa/(Fa*fe);
//     std::complex<double> alpha_0, a0_squared;
//     double r, a0_dagger_a0;
    
    
//     std::complex<double> M(double E, int beta); // Second order matrix element
// };
