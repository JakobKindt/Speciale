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
    // for (int i = 0; i < (int)taus.size(); ++i){
    //     std::cout << taus[i] << " ";
    // }
    // std::cout << "\n";
    // std::vector<double> taus{0, pi/10., pi*2/10., pi*3/10., pi*4/10., pi*5/10., pi*6/10., pi*7/10., pi*8/10., pi*9/10., pi*10/10., pi*11/10., pi*12/10., pi*13/10., pi*14/10., pi*15/10., pi*16/10., pi*17/10., pi*18/10., pi*19/10., pi*20/10.};
    double n = 1e16, varepsilon_0 = 0.5, r = std::asinh(std::sqrt(varepsilon_0*n));
    
    std::vector<double> angles{0., 0., pi/2, pi}, modula{0., r, r, r};

    // std::complex<double> zeta = std::polar(0., 0.);
    // sim S(zeta);
    // std::string outfile = "data.txt", outfile2 = "intensities.txt";
    // std::ofstream myoutput(outfile), myoutput2(outfile2);
    
    // for (int i = 0; i < (int)angles.size(); ++i){
    //     zeta = std::polar(modula[i], angles[i]); S.update(zeta);
    //         for (int j = 0; j < (int)Es.size(); ++j){
    //         E = Es[j];
    //         std::cout << "Test = " << S.test << ". Calculated a^2 = " << S.a0_squared << ", theoretical = " << S.alpha_0*S.alpha_0 - std::polar(1., std::arg(S.zeta))/2.*std::sinh(2*r)  << ". Calculated a^dagger a = " << S.a0_dagger_a0 << ", theoretical = " << std::norm(S.alpha_0) + std::sinh(S.r)*std::sinh(S.r) << "\n";
    //         myoutput << E << " " << S.A(E)*0.7729 << " " << S.C(E) << " " << S.theta(E) << "\n";
    //         }
    //         myoutput << "\n\n";
    //         for (int k = 0; k < (int)taus.size(); ++k){
    //             tau = taus[k];
    //             myoutput2 << tau;
    //             for (int j = 0; j < (int)E2s.size(); ++j){
    //             E = E2s[j];
                
    //             myoutput2 << " " << S.I(E, tau/omega_0);
    //             }
    //             myoutput2 << "\n";
    //         }
    //         myoutput2 << "\n\n";
    // }
    // myoutput.close();
    // myoutput2.close();


    std::complex<double> zeta = std::polar(r, 0.);
    // r = std::abs(z);
    std::complex<double> alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0), z = zeta;
    std::complex<double> phase = std::polar(1., std::arg(zeta));
    std::function<std::complex<double>(std::complex<double>)> Q = [&](std::complex<double> beta){
            std::complex<double> conj_beta = std::conj(beta), conj_alpha_0 = std::conj(alpha_0);
            return 1/(pi*std::cosh(r))*std::exp(-(std::norm(alpha_0) + std::norm(beta)) + (conj_beta*alpha_0 + beta*conj_alpha_0)/std::cosh(r) - ((phase*(conj_beta*conj_beta - conj_alpha_0*conj_alpha_0) + std::conj(phase)*(beta*beta - alpha_0*alpha_0))*std::tanh(r)/2.));
        };
    // std::cout << r << "\n";
    // Note that norm(z) returns norm squared of z, while abs(z) returns the norm.
    double L = (std::abs(alpha_0) + std::exp(r) + 1);
    std::function<std::complex<double>(std::complex<double>)> F = [=](std::complex<double> z){return std::norm(z)*std::exp(-std::norm(z));};
    std::complex<double> test = C_integrate_CC(Q, r - L, r + L);
    std::cout << "Testing integral = " << test << ", theoretical = " << pi << "\n";


    // auto test2 = [](double x) {
    //     return std::exp(-x*x);
    // };

    // std::cout << integrate_CC(test2, -10, 10, 1e-16, 1e-8) << '\n';

    // auto test3 = [](std::complex<double> z) {
    //     return std::exp(-std::norm(z));
    // };

    // std::cout << C_integrate_CC(test3, -inf, inf, 1e-8, 1e-8) << '\n';
return 0;
}