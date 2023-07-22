#include "mfem.hpp"
#include "cmath"

namespace mfem
{
    double t_final = 1.0;
    double dt = 0.001;

    //problem parameters
    double V = 1.0; //anode voltage
    double epsilon = 1.0;//permittivity

    ConstantCoefficient mu_e(1.0); //Mobility coefficents 20
    ConstantCoefficient mu_p(1.0);  // 0.0001

    ConstantCoefficient diff_const_e(1.0); //diffusion coefficients 200
    ConstantCoefficient diff_const_p(1.0); //.00005

    //both mobility and diffusion coefficients may not be constant depending on problem domain
    //should be possible to redefine variables as function coefficients 

    //steady state gas coefficients
    ConstantCoefficient Alpha(1); //ionization coefficient 200
    ConstantCoefficient Eta(0); //attachment coefficient 
    ConstantCoefficient Beta(1); //recombination coefficient; not used in perturbation equations
    ConstantCoefficient Gamma(0.1); //secondary electron emission coefficient 0.01

    ConstantCoefficient zero(0.0);
    ConstantCoefficient one(1.0);

}