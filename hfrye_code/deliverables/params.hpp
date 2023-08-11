#include "../../mfem.hpp"
#include "cmath"

namespace mfem
{
    const char *mesh_file = "../../data/inline-quad.mesh"; //Currently 2-d basic mesh to easily view results in visit
    //const char *mesh_file = "../../data/square-disc.mesh";
    //const char *mesh_file = "corner_plane.mesh"; //modified l-shape.mesh
    //const char *mesh_file = "corner_plane2.mesh";
    //const char *mesh_file = "../../data/amr-quad.mesh";
    //const char *mesh_file = "inline-segment.mesh";

    int order = 2; //second order legendre-Gauss solver (Make sure that's what it uses as according to Shibata paper)

    //Time evolution settings
    double t_final = 3.0;
    double dt = 0.01;
    const char *FE_db_name = "FE_evolution";

    bool run_fe_evolution = true;
    bool run_time_indep = false;

    bool direct_solve = false;

    //problem parameters
    double V = 37500; //anode voltage
    double epsilon = 8.8542e-12;//permittivity

    double pressure = 1013.2; //1013.2; //pressure mbar at atmospheric 

    //ConstantCoefficient mu_e(1.0); //Mobility coefficents 20
   // ConstantCoefficient mu_p(-1.0);  // 0.0001

   // ConstantCoefficient diff_const_e(1.0); //diffusion coefficients 200
   // ConstantCoefficient diff_const_p(1.0); //.00005

    //both mobility and diffusion coefficients may not be constant depending on problem domain
    //should be possible to redefine variables as function coefficients 

    //steady state gas coefficients; constant values
    //ConstantCoefficient Alpha(1.0); //ionization coefficient 200
    // ConstantCoefficient Eta(1.0); //attachment coefficient 
    //ConstantCoefficient Beta(0.0); //recombination coefficient; not used in perturbation equations
    
    ConstantCoefficient Gamma(0.1); //secondary electron emission coefficient 0.01 ; could potentially have its own
                                    //dependance on E_mag
    ConstantCoefficient zero(0.0);
    ConstantCoefficient one(1.0);

    //data fitted coefficient functions from BOLSIG+ data in Schnyder paper 

    double piecewise_func(double e){//for mu_e swarm parameter -- Schnyder
        
        if((e/pressure < 0.18) && (e/pressure > 0.0)){
            return 948.5 * tanh(-3.494 - 0.6097 * log(e / 100 * pressure)) + 948.5;
        } else if((e/pressure >= 0.18) && (e/pressure < 1150) ){
            return 139.6*pow(e/(100*pressure),-0.3455);
        } else if(e/pressure >= 1150){
            return 60;
        } else {
            return 0;
        }

    }

    //interpolated attachment coeeficient

    // interpolated wrt E field (linear)
    double interpolate_eta(double e){ //goes up to |E|=1500; here e is the horizontal value
        std::vector<double> eta_x = {0,20,40,60,100,150,200,300,400,500,700,1000,1200,1500};
        std::vector<double> eta_y = {0,0,.25,2.0,7.69,10,9.5,7.56,6.05,4.96,3.6,2.5,2.06,1.63};

        int size = eta_x.size();

        int i = 0;

        if(e >= eta_x[size-2]){
            i = size - 2;
        } else{
            while( e > eta_x[i+1]) i++;
        }
        double xL = eta_x[i], yL = eta_y[i], xR = eta_x[i+1], yR = eta_y[i+1];

        //extrapolation is used by default 
        if(e < xL) yR = yL;
        if(e > xR) yL = yR;

        double dydx = (yR - yL) / (xR - xL);

        return yL + dydx * ( e - xL );
        
    }

    
    //BCs for different meshes

    //returns the correct boundary markers for the cathode and anode of the 1D problem
    void lap_1d_bc(Array<int> &anode, Array<int> &cathode,Array<int> &ess_bdr){
        ess_bdr = 0;
        anode = 0;
        cathode = 0;

        anode[0] = 1;
        cathode[1] = 1;

        ess_bdr =1;

    }

    //returns the correct boundary markers for the cathode and anode of the 1D problem 
    //TODO: add Neumann conditions
    void lap_2d_bc(Array<int> &anode, Array<int> &cathode, Array<int> &ess_bdr){

        anode = 0;
        cathode = 0;
        ess_bdr = 0;

        anode[3] = 1;
        cathode[1] = 1;

        ess_bdr[3] = 1;
        ess_bdr[1] = 1;


    }

    void l_shape_bc(Array<int> &anode, Array<int> &cathode, Array<int> &ess_bdr){ //Have to modify boundary attr in example l-shape.mesh file first

        cathode = 0;
        anode = 0;

        cathode[0]=1;   
        cathode[1]=1;
        anode[3]=1;
        anode[4]=1;
        anode[5]=1; 

        ess_bdr = 0;

        ess_bdr[0]=1;
        ess_bdr[1]=1;
        ess_bdr[3]=1;
        ess_bdr[4]=1;
        ess_bdr[5]=1;

    }

    /*
    vertices
    10
    2
    0 0 
    0.2 0
    0.4 0
    0 0.05
    0.2 0.05
    0.4 0.05
    0.2 0.1
    0.4 0.1
    0.2 0.15
    0.4 0.15
    */

    void corner_plane_bc(Array<int> &anode, Array<int> &cathode, Array<int> &ess_bdr){
        cathode = 0;
        anode = 0;

        cathode[0]=1;
        anode[2]=1;

        ess_bdr = 0;

        ess_bdr[0] = 1;
        ess_bdr[2] = 1;
    }

    void lap_sqr_disc(Array<int> &anode, Array<int> &cathode, Array<int> &ess_bdr){

        cathode = 0 ;
        anode = 0;
        cathode[0] = 1;
        cathode[1] = 1;
        cathode[2] = 1;
        cathode[3] = 1;
        anode[4] = 1;
        anode[5] = 1;
        anode[6] = 1;
        anode[7] = 1;

        ess_bdr = 1;

    }


    //initial conditions for perturbations to FE evolution

    void velocity_function(const Vector &x, Vector &v)
    {
    int dim = x.Size();

    if(dim == 1){
    v(0) = 1.0;
    } else if(dim == 2){
        v(0) = 1.0;
        v(1) = 1.0;
    } else {
        std::cout << "dimension size not implmeneted yet. Aborting..."<<std::endl;
    }
    }

    double time_dep_diff_ic(const Vector &x){
        int dim = x.Size();
        Vector X(dim);

        if(dim == 1){
            return 1.0;
        } else if(dim == 2){
            //v(0) = sin(M_pI*X(0));
            //v(1) = sin(M_pI*X(1));
            //return  sin(M_pI*x(0)) + sin(M_pI*x(1));
            //return 0.0;
            return exp((-pow(x(0)-0.5,2)/0.1)-(pow(x(1)-0.5,2)/0.1));
        } else {
        std::cout << "dimension size not implmeneted yet. Aborting..."<<std::endl;
        return 1;
        }
    }

    double noise_ic(const Vector &x){
        int dim = x.Size();
        //Vector X(dim);
        
        unsigned int seed = time(NULL)+round(1000*x(0))+round(10000*x(1));
        srand(seed);

        double rand_num = ((double) rand() / (RAND_MAX));

        seed++;
        
        return rand_num;

    }


}