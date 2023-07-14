#include "mfem.hpp"
#include "cmath"

namespace mfem
{
    double t_final = 10.0;
    double dt = 0.01;

    //problem parameters
    double V = 1.0; //anode voltage
    double epsilon = 1.0;//permittivity

    ConstantCoefficient mu_e(1.0); //Mobility coefficents
    ConstantCoefficient mu_p(1.0); 

    ConstantCoefficient diff_const_e(1.0); //diffusion coefficients
    ConstantCoefficient diff_const_p(1.0);

    //both mobility and diffusion coefficients may not be constant depending on problem domain

    //steady state gas coefficients
    ConstantCoefficient Alpha(1.0); //ionization coefficient
    ConstantCoefficient Eta(1.0); //attachment coefficient
    ConstantCoefficient Beta(1.0); //recombination coefficient; not used in perturbation equations
    ConstantCoefficient Gamma(0.0); //secondary electron emission coefficient

    ConstantCoefficient zero(0.0);

}