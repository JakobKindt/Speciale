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
    
    // std::vector<double> angles{0., 0., pi/2, pi}, modula{0., r, r, r};

    // std::vector<double> angles{0}, modula{r};

    // std::complex<double> zeta = std::polar(0., 0.);
    // sim S(zeta);
    // std::string outfile = "data.txt", outfile2 = "intensities.txt";
    // std::ofstream myoutput(outfile), myoutput2(outfile2);
    
    // for (int i = 0; i < (int)angles.size(); ++i){
    //     zeta = std::polar(modula[i], angles[i]); S.update(zeta);
    //         for (int j = 0; j < (int)Es.size(); ++j){
    //         E = Es[j];
    //         // std::cout << "Calculated a^2 = " << S.a0_squared << ", theoretical = " << S.alpha_0*S.alpha_0 - std::polar(1., std::arg(S.zeta))/2.*std::sinh(2*r)  << ". Calculated a^dagger a = " << S.test << ", theoretical = " << std::norm(S.alpha_0) + std::sinh(S.r)*std::sinh(S.r) << "\n";
    //         std::cout << "Calculated a^2 = " << S.a0_squared << ", theoretical = " << S.alpha_0*S.alpha_0 - std::polar(1., std::arg(S.zeta))/2.*std::sinh(2*r)  << ". Calculated a^dagger a = " << S.a0_dagger_a0 << ", theoretical = " << std::norm(S.alpha_0) + std::sinh(S.r)*std::sinh(S.r) << "\n";
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


    // n = 1e1, varepsilon_0 = 0.5, r = std::asinh(std::sqrt(varepsilon_0*n));
    // double theta = 0;
    // std::complex<double> zeta = std::polar(r, theta);
    
    // std::complex<double> alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0), conj_alpha_0 = std::conj(alpha_0);
    // std::complex<double> alpha_0_squared = alpha_0*alpha_0, conj_alpha_0_squared = std::conj(alpha_0_squared);
    // std::complex<double> phase = std::polar(1., std::arg(zeta)), conj_phase = std::conj(phase);
    // double cosh_r = std::cosh(r), tanh_r = std::tanh(r), norm_alpha_0 = std::norm(alpha_0), prefactor = 1/(pi*cosh_r), inv_cosh_r = 1./cosh_r, half_tanh_r = tanh_r/2.;
    // double su = 1/std::sqrt(1 + tanh_r), sv = 1/std::sqrt(1 - tanh_r), susv = su*sv;
    // double x0 = std::abs(alpha_0) * std::exp(-r), L = 5, M = 3.25;
    // long int ncalls = 0;
    // auto Q = [&](std::complex<double> beta){
    //     ++ncalls;
    //     if (ncalls % (int)1e8 == 0){std::cerr << "ncalls = " << ncalls << "\n";}
    //     std::complex<double> conj_beta = std::conj(beta);
    //     return std::real(prefactor*std::exp(-(norm_alpha_0 + beta*conj_beta) + (conj_beta*alpha_0 + beta*conj_alpha_0)*inv_cosh_r - ((phase*(conj_beta*conj_beta - conj_alpha_0_squared) + conj_phase*(beta*beta - alpha_0_squared))*half_tanh_r)));
    // };


    // auto F = [&](std::complex<double> z){
    //     std::complex<double> beta = x0 + su*std::real(z) + std::complex<double>{0., 1.}* sv*std::imag(z);
    //     // return Q(beta)*susv*(std::norm(beta + alpha_0) - 1);
    //     return Q(beta)*susv;
    // };

    // std::complex<double> test = C_integrate_CC(F, x0 - M, x0 + M, -L, L, 1e2);
    // std::cout << "Testing integral = " << test << ", theoretical = " << 1 << ". After ncalls = " << ncalls << "\n";

    // std::complex<double> gamma{1.2, 0.7};
    // std::cout << "F(gamma) = " << F(gamma) << "\n";


    
    // std::complex<double> zeta = std::polar(r, 0.);
    // double theta = std::arg(zeta);
    // std::complex<double> alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
    // double xa = std::real(alpha_0), ya = std::imag(alpha_0);
    // std::complex<double> alpha_0_squared = alpha_0*alpha_0;
    // std::complex<double> phase = std::polar(1., theta);

    // double cosh_r = std::cosh(r), tanh_r = std::tanh(r), norm_alpha_0 = std::norm(alpha_0), prefactor = 1/(pi*cosh_r), inv_cosh_r_2 = 2./cosh_r, sin_tanh = std::sin(theta)*tanh_r, cos_tanh = std::cos(theta)*tanh_r;
    // double su = 1/std::sqrt(1 + tanh_r), sv = 1/std::sqrt(1 - tanh_r), susv = su*sv;
    // double exponent = std::real(phase*alpha_0_squared)*tanh_r;
    // // double exponent = - std::norm(alpha_0)*(1.0 - tanh_r) - (1.0 + tanh_r)*xa*xa - (1.0 - tanh_r)*ya*ya + std::abs(alpha_0)*xa*inv_cosh_r_2;
    // double x0 = std::abs(alpha_0) * std::exp(-r), L = 5, M = 3.25;
    // long int ncalls = 0;
    

    // auto Q = [&](std::complex<double> beta){
    //     ++ncalls;
    //     if (ncalls % (int)1e8 == 0){std::cerr << "ncalls = " << ncalls << "\n";}
    //     double x = std::real(beta), y = std::imag(beta);
    //     return prefactor*std::exp(-(norm_alpha_0 + x*x + y*y) + (x*xa + y*ya)*inv_cosh_r_2 + exponent + 2*sin_tanh*x*y + cos_tanh*(y*y - x*x));
    // };


    // auto F = [&](std::complex<double> z){
    //     std::complex<double> beta = x0 + su*std::real(z) + std::complex<double>{0., 1.}* sv*std::imag(z);
    //     return Q(beta)*susv;
    // };

    // std::complex<double> test = C_integrate_CC(F, x0 - M, x0 + M, -L, L, 1e-5);
    // std::cout << "Testing integral = " << test << ", theoretical = " << 1 << ". After ncalls = " << ncalls << "\n";

    // // std::complex<double> gamma{1.2, 0.7};
    // // std::cout << "F(gamma) = " << F(gamma) << "\n";






    // std::complex<double> zeta = std::polar(r, 0.);
    // double theta = std::arg(zeta);
    // std::complex<double> alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
    // double cosh_r = std::cosh(r), tanh_r = std::tanh(r), prefactor = 1/(pi*cosh_r), sin_tanh = std::sin(theta)*tanh_r, cos_tanh = std::cos(theta)*tanh_r;
    // double su = 1/std::sqrt(1 + tanh_r), sv = 1/std::sqrt(1 - tanh_r), susv = su*sv;
    // double x0 = std::abs(alpha_0) * std::exp(-r), L = 5, M = 3.25;
    // long int ncalls = 0;
    

    // auto Q = [&](std::complex<double> beta){
    //     ++ncalls;
    //     if (ncalls % (int)1e8 == 0){std::cerr << "ncalls = " << ncalls << "\n";}
    //     double x = std::real(beta), y = std::imag(beta);
    //     // return prefactor*std::exp(-(norm_alpha_0 + x*x + y*y) + (x*xa + y*ya)*inv_cosh_r_2 + exponent + 2*sin_tanh*x*y + cos_tanh*(y*y - x*x));
    //     return prefactor*std::exp(-(x*x + y*y) + 2*sin_tanh*x*y + cos_tanh*(y*y - x*x));
    // };


    // auto F = [&](std::complex<double> z){
    //     std::complex<double> beta = x0 + su*std::real(z) + std::complex<double>{0., 1.}* sv*std::imag(z);
    //     // return Q(beta)*susv*(std::norm(beta + alpha_0) - 1);
    //     return Q(beta)
    // };

    // std::complex<double> test = C_integrate_CC(F, x0 - M, x0 + M, -L, L, 1e2);
    // std::cout << "Testing integral = " << test << ", theoretical = " << 1 << ". After ncalls = " << ncalls << "\n";

    // std::complex<double> gamma{1.2, 0.7};
    // std::cout << "F(gamma) = " << F(gamma) << "\n";








    // std::complex<double> zeta = std::polar(r, 0.);
    // double theta = std::arg(zeta);
    // std::complex<double> alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
    // double cosh_r = std::cosh(r), tanh_r = std::tanh(r), prefactor = 1/(pi*cosh_r), sin_tanh = std::sin(theta)*tanh_r, cos_tanh = std::cos(theta)*tanh_r;
    // double su = 1/std::sqrt(1 + tanh_r), sv = 1/std::sqrt(1 - tanh_r), susv = su*sv;
    // double x0 = std::abs(alpha_0) * std::exp(-r), L = 5, M = 3.25;
    // long int ncalls = 0;
    

    // auto Q = [&](std::complex<double> beta){
    //     ++ncalls;
    //     if (ncalls % (int)1e8 == 0){std::cerr << "ncalls = " << ncalls << "\n";}
    //     // double x = std::real(beta), y = std::imag(beta);
    //     // return prefactor*std::exp(-(norm_alpha_0 + x*x + y*y) + (x*xa + y*ya)*inv_cosh_r_2 + exponent + 2*sin_tanh*x*y + cos_tanh*(y*y - x*x));
    //     return std::exp(-std::norm(beta))/pi;
    // };


    // auto F = [&](std::complex<double> z){
    //     std::complex<double> beta = x0 + su*std::real(z) + std::complex<double>{0., 1.}* sv*std::imag(z);
    //     // return Q(beta)*susv*(std::norm(beta) - 1);
    //     return Q(beta)*susv;
    // };

    // std::complex<double> test = C_integrate_CC(F, -M, M, -L, L, 1e2);
    // std::cout << "Testing integral = " << test << ", theoretical = " << 1 << ". After ncalls = " << ncalls << "\n";

    // std::complex<double> gamma{1.2, 0.7};
    // std::cout << "F(gamma) = " << F(gamma) << "\n";










// auto Q = [&](std::complex<double> beta){
    //     ++calls;
    //     double x = std::real(Q), y = std::imag(Q);
    //     std::complex<double> conj_beta = std::conj(beta);
    //     return std::real(prefactor*std::exp(-(norm_alpha_0 + beta*conj_beta) + (conj_beta*alpha_0 + beta*conj_alpha_0)*inv_cosh_r - ((phase*(conj_beta*conj_beta - conj_alpha_0_squared) + conj_phase*(beta*beta - alpha_0_squared))*half_tanh_r)));
    // };




    double theta = pi/2.;
    std::complex<double> zeta = std::polar(r, theta);
    
    std::complex<double> alpha_0 = std::polar(std::sqrt(n*(1.0-varepsilon_0)), 0.0);
    // double re_alpha_0 = std::real(alpha_0), im_alpha_0 = std::imag(alpha_0);
    std::complex<double> phase = std::polar(1., std::arg(zeta));
    double cosh_r = std::cosh(r), tanh_r = std::tanh(r);
    double su = 1/std::sqrt(1 + tanh_r), sv = 1/std::sqrt(1 - tanh_r); // susv = su*sv = cosh(r), so prefactor*cosh(r) = 1/pi which is more nummerical stable;
    su = std::sqrt(0.5*(std::exp(-2.0*r) + 1.0)); sv = std::sqrt(0.5*(std::exp(2.0*r) + 1.0)); // ChatGPT. Verify later
    // double x0 = std::abs(alpha_0) * std::exp(-r), L = 5, M = 3.25;
    // double L = 5, M = 3.25;
    double L = 10, M = 10;
    long int ncalls = 0;
    auto Q = [&](std::complex<double> beta){
        ++ncalls;
        if (ncalls % (int)1e8 == 0){std::cerr << "ncalls = " << ncalls << "\n";}
        // std::complex<double> conj_beta = std::conj(beta);
        return std::exp(-std::norm(beta))/pi;
    };
    
    std::complex<double> re, im;
    double cos_2 = std::cos(theta/2), sin_2 = std::sin(theta/2);

    auto F_norm = [&](std::complex<double> z){
        re = std::real(z); im = std::imag(z);
        std::complex<double> beta = alpha_0 + su*(cos_2*re + sin_2*im) + std::complex<double>{0., 1.}* sv*(-sin_2*re + cos_2*im);
        return Q(z);
    };

    auto F = [&](std::complex<double> z){
        // std::complex<double> beta = su*(std::real(z) + re_alpha_0) + std::complex<double>{0., 1.}* sv*(std::imag(z) + im_alpha_0);
        // std::complex<double> beta = alpha_0 + su*std::real(z) + std::complex<double>{0., 1.}* sv*std::imag(z);
        re = std::real(z); im = std::imag(z);
        std::complex<double> beta = alpha_0 + su*(cos_2*re + sin_2*im) + std::complex<double>{0., 1.}* sv*(-sin_2*re + cos_2*im);
        return Q(z)*(std::norm(beta) - 1);
        // return Q(z)*susv;
    };
    
    
    auto F2 = [&](std::complex<double> z){
        // std::complex<double> beta = su*(std::real(z) + re_alpha_0) + std::complex<double>{0., 1.}* sv*(std::imag(z) + im_alpha_0);
        re = std::real(z); im = std::imag(z);
        std::complex<double> beta = alpha_0 + su*(cos_2*re + sin_2*im) + std::complex<double>{0., 1.}* sv*(-sin_2*re + cos_2*im);
        return Q(z)*beta*beta;
        // return Q(z)*std::real(beta)*std::real(beta);
        // return Q(z)*std::imag(beta)*std::imag(beta);
    };
    int N = 1e3;
    std::complex<double> test_norm = C_integrate_CC(F_norm, -M,  M, -L, L, N, N);
    std::complex<double> test = C_integrate_CC(F, -M,  M, -L, L, N, N);
    std::cout << "Testing norm = " << test_norm << ", theoretical norm = 1. After ncalls = " << ncalls << "\n";
    std::cout << "Testing integral = " << test << ", theoretical a_dagger a = " << std::norm(alpha_0) + std::sinh(r)*std::sinh(r) << ". After ncalls = " << ncalls << "\n";
    std::complex<double> test2 = C_integrate_CC(F2, -M,  M, -L, L, N*1e2, N*1e2);
    std::cout << "Testing integral = " << test2 << ", theoretical a^2 = " << alpha_0*alpha_0 - phase*std::sinh(2*r)/2. << ". After ncalls = " << ncalls << "\n";
    std::cout << "Theoretical a = " << alpha_0 << ".\n";


    std::cout << "su = " << su << ". sv = " << sv << "\n";


    
return 0;
}


// std::cout << r << "\n";
// Note that norm(z) returns norm squared of z, while abs(z) returns the norm.
// double L = 2*(std::abs(alpha_0) + std::exp(r) + 1);
// std::function<std::complex<double>(std::complex<double>)> F = [=](std::complex<double> z){return std::norm(z)*std::exp(-std::norm(z));};
// L = 5e8;
// double x0 = std::abs(alpha_0) * std::exp(-r);
// double M = 1.8;