#include"Bunto.h"
// #include"hf.h"
#include<fstream>
#include<iostream>

// int main(int argc, char** argv)
int main(){
    std::cerr << "Begun \n";
return 0;
}


// #include"Bunto.h"
// #include"hf.h"
// #include<fstream>
// #include<iostream>

// // int main(int argc, char** argv)
// int main(){
//     std::cerr << "Begun \n";
//     std::vector<double> Es{1.892, 4.990, 8.088, 11.19, 14.29, 17.38, 20.48, 23.58};
//     std::vector<double> E2s{1.892, 8.088, 14.29, 20.48};
//     double E, tau, pi = 3.14159265358979323846, omega_0 = 0.057;
//     std::vector<double> taus = linspace(0., 2*pi, 40);
//     double n = 1e16, varepsilon_0 = 0.5, r = std::asinh(std::sqrt(varepsilon_0*n));
    
//     // std::vector<double> angles{0., 0., pi/2, pi}, modula{0., r, r, r};

//     // std::vector<double> angles{0., pi}, modula{0., r};
//     std::vector<double> angles{0}, modula{r};

//     // std::complex<double> zeta = std::polar(0., 0.);
//     // sim S(zeta);
//     // std::string outfile = "data.txt", outfile2 = "intensities.txt";
//     // std::ofstream myoutput(outfile), myoutput2(outfile2);
    
//     // for (int i = 0; i < (int)angles.size(); ++i){
//     //     zeta = std::polar(modula[i], angles[i]); S.update(zeta);
//     //         for (int j = 0; j < (int)Es.size(); ++j){
//     //         E = Es[j];
//     //         // std::cout << "Calculated a^2 = " << S.a0_squared << ", theoretical = " << S.alpha_0*S.alpha_0 - std::polar(1., std::arg(S.zeta))/2.*std::sinh(2*r)  << ". Calculated a^dagger a = " << S.test << ", theoretical = " << std::norm(S.alpha_0) + std::sinh(S.r)*std::sinh(S.r) << "\n";
//     //         std::cout << "Calculated a^2 = " << S.a0_squared << ", theoretical = " << S.alpha_0*S.alpha_0 - std::polar(1., std::arg(S.zeta))/2.*std::sinh(2*r)  << ". Calculated a^dagger a = " << S.a0_dagger_a0 << ", theoretical = " << std::norm(S.alpha_0) + std::sinh(S.r)*std::sinh(S.r) << "\n";
//     //         myoutput << E << " " << S.A(E)*0.7729 << " " << S.C(E) << " " << S.theta(E) << "\n";
//     //         }
//     //         myoutput << "\n\n";
//     //         for (int k = 0; k < (int)taus.size(); ++k){
//     //             tau = taus[k];
//     //             myoutput2 << tau;
//     //             for (int j = 0; j < (int)E2s.size(); ++j){
//     //             E = E2s[j];
                
//     //             myoutput2 << " " << S.I(E, tau/omega_0);
//     //             }
//     //             myoutput2 << "\n";
//     //         }
//     //         myoutput2 << "\n\n";
//     // }
//     // myoutput.close();
//     // myoutput2.close();




//     std::complex<double> zeta = std::polar(r, 0.);
    
//     // r = std::abs(z);
//     std::complex<double> alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
//     std::complex<double> phase = std::polar(1., std::arg(zeta));

//     double su = 1/std::sqrt(1 + std::tanh(r)), sv = 1/std::sqrt(1 - std::tanh(r));
//     double x0 = std::abs(alpha_0) * std::exp(-r), L = 5, M = 3.25;
//     int calls = 0;
//     std::cerr << "Before \n";
//     auto Q = [&](std::complex<double> beta){
//         std::complex<double> conj_beta = std::conj(beta), conj_alpha_0 = std::conj(alpha_0);
//         ++calls;
//         std::cerr << calls << "\n";
//         return 1/(pi*std::cosh(r))*std::exp(-(std::norm(alpha_0) + std::norm(beta)) + (conj_beta*alpha_0 + beta*conj_alpha_0)/std::cosh(r) - ((phase*(conj_beta*conj_beta - conj_alpha_0*conj_alpha_0) + std::conj(phase)*(beta*beta - alpha_0*alpha_0))*std::tanh(r)/2.));
//     };

//     auto F = [&](std::complex<double> z){
//         std::complex<double> beta = x0 + su*std::real(z) + std::complex<double>{0., 1.}* sv*std::imag(z);
//         return Q(beta)*su*sv;
//     };

    
//     // std::cout << r << "\n";
//     // Note that norm(z) returns norm squared of z, while abs(z) returns the norm.
//     // double L = 2*(std::abs(alpha_0) + std::exp(r) + 1);
//     // std::function<std::complex<double>(std::complex<double>)> F = [=](std::complex<double> z){return std::norm(z)*std::exp(-std::norm(z));};
//     // L = 5e8;
//     // double x0 = std::abs(alpha_0) * std::exp(-r);
//     // double M = 1.8;
//     std::complex<double> test = C_integrate_CC(F, x0 - M, x0 + M, -L, L, 1e30, 1e30);
//     // std::complex<double> test = C_integrate_CC(F, x0 - M, x0 + M, -L, L, 1e-5);
//     std::cout << "Testing integral = " << test << ", theoretical = " << 1 << "\n";
// return 0;
// }