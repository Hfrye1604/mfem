#include "../../mfem.hpp"
#include "cmath"

namespace mfem
{
    double t_final = 1.0;
    double dt = 0.001;

    //problem parameters
    double V =1.0; //anode voltage
    double epsilon = 1.0;//permittivity

    ConstantCoefficient mu_e(1.0); //Mobility coefficents 20
    ConstantCoefficient mu_p(-1.0);  // 0.0001

    ConstantCoefficient diff_const_e(1.0); //diffusion coefficients 200
    ConstantCoefficient diff_const_p(1.0); //.00005

    //both mobility and diffusion coefficients may not be constant depending on problem domain
    //should be possible to redefine variables as function coefficients 

    //steady state gas coefficients
    ConstantCoefficient Alpha(.7); //ionization coefficient 200
    ConstantCoefficient Eta(.3); //attachment coefficient 
    ConstantCoefficient Beta(1); //recombination coefficient; not used in perturbation equations
    ConstantCoefficient Gamma(10.0); //secondary electron emission coefficient 0.01

    ConstantCoefficient zero(0.0);
    ConstantCoefficient one(1.0);

    //data fitted coefficient functions from BOLSIG+ data in Schnyder paper alpha/p
    //Expressions are computed divided by pressure */p
/*
    double Alpha_func(double V_applied){
        double A = 2900;
        double B = 28070;//expressions uses from table 14.1 in Lieberman paper found through schynder paper
        double p = 1.0; 

        return A * exp((-B*p)/V_applied);
    }

    double Eta_func(double V_applied){
        double A = .0000009;
        double B = 200;

        return 
    }
    */

}