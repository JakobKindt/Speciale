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
    std::complex<double> zeta;
    sim(std::complex<double> z){
        zeta = z; 
        r = std::abs(z);
        alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
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
    std::complex<double> test;
    
    std::complex<double> M(double E, int beta); // Second order matrix element
};