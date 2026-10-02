#include"Bunto.h"
#include"hf.h"
#include<fstream>
#include<iostream>

// int main(int argc, char** argv)
int main(){
    std::vector<double> Es{1.892, 4.990, 8.088, 11.19, 14.29, 17.38, 20.48, 23.58};
    std::vector<double> E2s{1.892, 8.088, 14.29, 20.48};
    double E, tau, pi = 3.14159265358979323846, omega_0 = 0.057;
    std::vector<double> taus = linspace(0., 2*pi, 40);
    double n = 1e16, varepsilon_0 = 0.5, r = std::asinh(std::sqrt(varepsilon_0*n));
    
    std::vector<double> angles{0., 0., pi/2, pi}, modula{0., r, r, r};

    std::complex<double> zeta = std::polar(0., 0.);
    sim S(zeta);
    std::string outfile = "data.txt", outfile2 = "intensities.txt";
    std::ofstream myoutput(outfile), myoutput2(outfile2);
    
    for (int i = 0; i < (int)angles.size(); ++i){
        zeta = std::polar(modula[i], angles[i]); S.update(zeta);
            for (int j = 0; j < (int)Es.size(); ++j){
            E = Es[j];
            // std::cout << "Calculated a^2 = " << S.a0_squared << ", theoretical = " << S.alpha_0*S.alpha_0 - std::polar(1., std::arg(S.zeta))/2.*std::sinh(2*r)  << ". Calculated a^dagger a = " << S.test << ", theoretical = " << std::norm(S.alpha_0) + std::sinh(S.r)*std::sinh(S.r) << "\n";
            myoutput << E << " " << S.A(E)*0.7729 << " " << S.C(E) << " " << S.theta(E) << "\n";
            }
            myoutput << "\n\n";
            for (int k = 0; k < (int)taus.size(); ++k){
                tau = taus[k];
                myoutput2 << tau;
                for (int j = 0; j < (int)E2s.size(); ++j){
                E = E2s[j];
                
                myoutput2 << " " << S.I(E, tau/omega_0);
                }
                myoutput2 << "\n";
            }
            myoutput2 << "\n\n";
        std::cout << "Calculated a^2 = " << S.a0_squared << ", theoretical = " << S.alpha_0*S.alpha_0 - std::polar(1., std::arg(S.zeta))/2.*std::sinh(2*r)  << ". Calculated a^dagger a = " << S.a0_dagger_a0 << ", theoretical = " << std::norm(S.alpha_0) + std::sinh(S.r)*std::sinh(S.r) << "\n";
    }
    myoutput.close();
    myoutput2.close();





    // double theta = 0;
    // std::complex<double> zeta = std::polar(r, theta);
    
    // std::complex<double> alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
    // std::complex<double> phase = std::polar(1., std::arg(zeta));
    // std::complex<double> phase_2 = std::polar(1., std::arg(zeta)/2);
    // // std::complex<double> phase{-1., 0.};
    // double tanh_r = std::tanh(r), su = 1/std::sqrt(1 + tanh_r), sv = 1/std::sqrt(1 - tanh_r); // susv = su*sv = cosh(r), so prefactor*cosh(r) = 1/pi which is more nummerical stable;
    // su = std::sqrt(0.5*(std::exp(-2.0*r) + 1.0)); sv = std::sqrt(0.5*(std::exp(2.0*r) + 1.0)); // ChatGPT. Verify later
    // // double x0 = std::abs(alpha_0) * std::exp(-r), L = 5, M = 3.25;
    // // double L = 5, M = 3.25;
    // double L = 10, M = 10;
    // long int ncalls = 0;
    // auto Q = [&](std::complex<double> beta){
    //     ++ncalls;
    //     if (ncalls % (int)1e8 == 0){std::cerr << "ncalls = " << ncalls << "\n";}
    //     return std::exp(-std::norm(beta))/pi;
    // };
    
    // std::complex<double> re, im;
    // // double cos_2 = std::cos(theta/2.), sin_2 = std::sin(theta/2.);
    // double cos_2 = std::cos(theta/2.), sin_2 = std::sin(theta/2.);
    // // double cos_2 = 0, sin_2 = 1;
    // std::cout << "cos_2 = " << cos_2 << ". sin_2 = " << sin_2 << ".\n";

    // auto F_norm = [&](std::complex<double> z){
    //     re = std::real(z); im = std::imag(z);
    //     return Q(z);
    // };

    // auto F = [&](std::complex<double> z){
    //     re = std::real(z); im = std::imag(z);
    //     std::complex<double> beta = alpha_0 + (su*cos_2*re - sv*sin_2*im) + std::complex<double>{0., 1.}*(su*sin_2*re + sv*cos_2*im);
    //     return Q(z)*(std::norm(beta) - 1);
    // };
    
    
    // auto F2 = [&](std::complex<double> z){
    //     re = std::real(z); im = std::imag(z);
    //     std::complex<double> beta = alpha_0 + (su*cos_2*re - sv*sin_2*im) + std::complex<double>{0., 1.}*(su*sin_2*re + sv*cos_2*im);
    //     return Q(z)*beta*beta;
    // };
    // int N = 1e3;
    // std::complex<double> test_norm = C_integrate_CC(F_norm, -M,  M, -L, L, N, N);
    // std::complex<double> test = C_integrate_CC(F, -M,  M, -L, L, N, N);
    // std::cout << "Measured norm = " << test_norm << ", theoretical norm = 1. After ncalls = " << ncalls << "\n";
    // std::cout << "Measured a_dagger a = " << test << ", theoretical a_dagger a = " << std::norm(alpha_0) + std::sinh(r)*std::sinh(r) << ". After ncalls = " << ncalls << "\n";
    // std::complex<double> test2 = C_integrate_CC(F2, -M,  M, -L, L, N*1e1, N*1e1);
    // std::cout << "Phase = " << phase << ".\n";
    // std::cout << "Measured a^2 = " << test2 << ", theoretical a^2 = " << alpha_0*alpha_0 - phase*std::sinh(2*r)/2. << ". After ncalls = " << ncalls << "\n";
    // std::cout << "Theoretical a = " << alpha_0 << ".\n";


    // std::cout << "su = " << su << ". sv = " << sv << "\n";

    
return 0;
}