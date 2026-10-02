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

    std::complex<double> zeta = std::polar(0., 0.);
    sim S(zeta);
    std::string outfile = "data.txt", outfile2 = "intensities.txt";
    std::ofstream myoutput(outfile), myoutput2(outfile2);
    
    for (int i = 0; i < (int)angles.size(); ++i){
        zeta = std::polar(modula[i], angles[i]); S.update(zeta);
            for (int j = 0; j < (int)Es.size(); ++j){
            E = Es[j];
            std::cout << "Calculated a^2 = " << S.a0_squared << ", theoretical = " << S.alpha_0*S.alpha_0 - std::polar(1., std::arg(S.zeta))/2.*std::sinh(2*r)  << ". Calculated a^dagger a = " << S.a0_dagger_a0 << ", theoretical = " << std::norm(S.alpha_0) + std::sinh(S.r)*std::sinh(S.r) << "\n";
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
    }
    myoutput.close();
    myoutput2.close();
    // std::cout << r << "\n";
    // Note that norm(z) returns norm squared of z, while abs(z) returns the norm.
return 0;
}