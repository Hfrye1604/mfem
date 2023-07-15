#include "mfem.hpp"
#include "cmath"

namespace mfem
{
    double t_final = 0.1;
    double dt = 0.001;

    //problem parameters
    double V = 0.00001; //anode voltage
    double epsilon = 1.0;//permittivity

    ConstantCoefficient mu_e(1.00); //Mobility coefficents
    ConstantCoefficient mu_p(1.00); 

    ConstantCoefficient diff_const_e(10); //diffusion coefficients
    ConstantCoefficient diff_const_p(10);

    //both mobility and diffusion coefficients may not be constant depending on problem domain
    //should be possible to redefine variables as function coefficients 

    //steady state gas coefficients
    ConstantCoefficient Alpha(0.5); //ionization coefficient
    ConstantCoefficient Eta(0.5); //attachment coefficient
    ConstantCoefficient Beta(1.0); //recombination coefficient; not used in perturbation equations
    ConstantCoefficient Gamma(0.1); //secondary electron emission coefficient

    ConstantCoefficient zero(0.0);
    ConstantCoefficient one(0.0);

}