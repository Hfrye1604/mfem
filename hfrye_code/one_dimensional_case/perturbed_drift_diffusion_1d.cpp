/*1-d plate to plate conduction discharge tube code using MFEM for the Drift-Diffusion equations
    --The purpose of this code is to attempt to solve the D-D perturbation equations 
    to verify the physics before defining for the jacobian and mass matrices of the system
    for stability analysis. This code is written with dimensional generality in mind.
*/
//#include "/g/g11/frye11/mfem/mfem.hpp"
#include "mfem.hpp"
//#include "NSnonlininteg_modified.hpp"
//#include "fem/advectionSUPGinteg.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

//Compute E field
Vector E_field(int, double, double, FiniteElementSpace);

//return correct boundary markers for problem
void lap_1d_bc(Array<int> &anode, Array<int> &cathode);
void lap_2d_bc(Array<int> &anode, Array<int> &cathode);

//velocity coefficient -- input to program based off of electric field E
void velocity_function_e(const Vector &x, Vector &v);
void velocity_function_p(const Vector &x, Vector &v);

//BC for secondary electron emission
double sec_emi_funct(const Vector &x);

//init conidtion -- non-ionized state
double ne0_funct(const Vector &x);
double np0_funct(const Vector &x);

//Boundary condition (might need to handle differently than in example 9)
void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();

   v(0) = 1.0;
}


int main(int agrc, char *argv[]){
    //const char *mesh_file = "../../mfem/data/inline-quad.mesh"; //Currently 2-d basic mesh to easily view results in visit
    //const char *mesh_file = "../../mfem/data/inline-segment.mesh";
    int order = 2; //second order legendre-Gauss solver (Make sure that's what it uses as according to Shibata paper)
    double t_final = 10.0;
    double dt = 0.01;

    //time integrator for model verification purposes, perturbation equations might carry little physical meaning integrated like this however.
    ODESolver *ode_solver = NULL;
    ode_solver = new BackwardEulerSolver;

    //Define Mesh 
    //Mesh mesh(mesh_file);

    Mesh mesh(10,1.0);
    // possibly refine mesh to increase resolution
    mesh.UniformRefinement();
   // mesh.UniformRefinement();

    int dim = mesh.Dimension();

    //define finite element space; for now, try to stick to H1 space and implement nonlinNSint but might need different space
    H1_FECollection fec(order,dim); //Gauss-Legendre unsupported
    FiniteElementSpace fespace(&mesh,&fec);
    cout << "Number of unknowns: " << fespace.GetTrueVSize() << endl;

    //boundary markers for cathode and anode
    //For 1D case
    Array<int> cathode_bdr(mesh.bdr_attributes.Max());
    Array<int> anode_bdr(mesh.bdr_attributes.Max());
    cathode_bdr =0;
    anode_bdr = 0;
    cathode_bdr[0]=1;
    anode_bdr[1]=1;
    //E_field for 1d

    double V = 1.0; //anode voltage
    Vector E; 
    E = E_field(order,V,1.0,fespace);
/* for diagnostic purposes
    DataCollection *dc = NULL;
    dc = new VisItDataCollection("E_field",&mesh);
    dc->SetPrecision(8);
    dc->RegisterField("solution",&E);
    dc->Save();
*/
    //apply BCs here
    //Array<int>

    //set up linear and bilinear forms for each 
    //set up velocity coefficients for particle species
    //VectorFunctionCoefficient velocity_e(dim, velocity_function_e);
    //VectorFunctionCoefficient velocity_p(dim, velocity_function_p);
    double mu_e = 1.0;//will add accurate values later
    double mu_p = 1.0;

    //GridFunctionCoefficient velocity_e(mu_e*E);
    //GridFunctionCoefficient velocity_p(mu_p*E);
    Vector v_e;
    v_e = E;
    v_e *= mu_e;

    cout << "v_e:" << endl;
    for(int i = 0; i < v_e.Size();i++)
        cout << v_e[i] << endl;

    VectorConstantCoefficient velocity_e(v_e);

    Vector v_p;
    v_p = E;
    v_p *= mu_p;

    Vector v_test;
    v_test.SetSize(1);
    v_test[0]=1.0;

    //Vector v_test(-0.1,1);

    VectorConstantCoefficient velocity_p(v_p);
    VectorFunctionCoefficient vel_test(dim,velocity_function);
    //cout << "vel_test" << vel_test.Size() << endl;
    //for(int i=0;i < vel_test.Size();i++)
      //  cout << vel_test[i] << endl;
        
    //BCs will have to specify in function the specifics of BCs later
    
    //For initial conditions primarily for physics check, will have to specify later

    ConstantCoefficient zero(0.0);
    
    //Following integrator matrix notation from Shibata paper
    BilinearForm M(&fespace); //mass matrix the same for both species
    BilinearForm Ae(&fespace);
    BilinearForm De(&fespace);
    BilinearForm Kee(&fespace);
    BilinearForm Kep(&fespace);
    BilinearForm See(&fespace);

    BilinearForm Ap(&fespace);
    BilinearForm Dp(&fespace);
    BilinearForm Kp(&fespace);
    BilinearForm Spe(&fespace);

    Vector vec;
    vec.SetSize(1);
    vec[0] = 1.0;

    VectorConstantCoefficient vtest(vec);

    ConstantCoefficient diff_const_e(1.0);
    ConstantCoefficient diff_const_p(1.0);

    M.AddDomainIntegrator(new MassIntegrator); //same for both equations
    
    ConstantCoefficient viscocity(0.0); 
    Ae.AddDomainIntegrator(new ConservativeConvectionIntegrator(vel_test,-1.0)); //alpha=-1.0 to flip sign
    Ae.AddDomainIntegrator(new AdvectionSUPGIntegrator(vtest,0.0,1.0));

    De.AddDomainIntegrator(new DiffusionIntegrator(diff_const_e));

    ConstantCoefficient gamma(0.1); //to be multipled into the mass matrix of the Kep operator
    Kee.AddBoundaryIntegrator(new MassIntegrator,cathode_bdr); //TODO: will have to verify correct boundary term ;will need to split up into Kee and Kpe terms
    //Kee.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(velocity_e), cathode_bdr);

    Kep.AddBoundaryIntegrator(new MassIntegrator(gamma),cathode_bdr);
    //Kep.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(velocity_p), cathode_bdr);

    ConstantCoefficient SeeConst(1.0); //TODO: will have to compute full expression constant 
    See.AddDomainIntegrator(new MassIntegrator(SeeConst));

    Ap.AddDomainIntegrator(new ConservativeConvectionIntegrator(vel_test,-1.0));
    Ap.AddDomainIntegrator(new AdvectionSUPGIntegrator(vtest,0.0,2.0));

    Dp.AddDomainIntegrator(new DiffusionIntegrator(diff_const_p));

    Kep.AddBoundaryIntegrator(new MassIntegrator,anode_bdr);
    //Kep.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(velocity_p), anode_bdr);

    ConstantCoefficient SpeConst(1.0); //TODO: will have to compute full See expression constant
    Spe.AddDomainIntegrator(new MassIntegrator(SpeConst));

    //after integrators are correctly defined, the Jacobian can be simply added or subtracted together 
    // as specified after assembling
    
    M.Assemble();
    Ae.Assemble();   
    De.Assemble();   
    //Kee.Assemble();
    //Kep.Assemble();
    See.Assemble();

    Ap.Assemble();
    Dp.Assemble();
    //Kp.Assemble();
    Spe.Assemble();

    // set up Jacobians

    SparseMatrix J11, AeSp, DeSp;
    AeSp = Ae.SpMat();
    //AeSp *= -1.0;
    DeSp = De.SpMat();
    //DeSp *= -1.0; 
    //J11 = -Ae.SpMat() - De.SpMat() + Kee.SpMat() + See.SpMat();
    J11.Add(-1.0,Ae.SpMat());
    J11.Add(-1.0,De.SpMat());
    J11.Add(1.0,See.SpMat());

    //cout << J11 << endl;

    SparseMatrix J12; 
   // J12 = Kep.SpMat();

    SparseMatrix J21; 
    

    SparseMatrix J22;
    //J22 = -Ap.SpMat() - Dp.SpMat() + Kp.SpMat();

    //Set up block jacobian and block mass matrix


    //eigensolver TBD

   //ofstream meshfile;
    //meshfile.open("mesh.dat");
    //for(int i=0; i < mesh.GetElementSize)

    return 0;
}

//solves preliminary laplace problem for the electric potential and then returns the resulting electric field
Vector E_field(int order, double V, double epsilon, FiniteElementSpace fes){// pass in boundary elements along with fes since it should match up with d-d problem 
//permittivity
    Mesh *mesh = fes.GetMesh();

    //Array<int> ess_bdr(mesh.bdr_attributes.Max());
     // ess_bdr = 1;
      //fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    //u represents the electric potential
    GridFunction u(&fes);
    u = 0.0;
    //u[0]=V;
    
    cout.precision(8);
    cout << "initialized :" << endl;
    cout << u << endl;
    
    ConstantCoefficient VCoef(V); //anode voltage
    ConstantCoefficient epCoef(epsilon);

    BilinearForm a(&fes);
    a.AddDomainIntegrator(new DiffusionIntegrator(epCoef));
    a.Assemble();

    ConstantCoefficient zero(0.0);
    LinearForm b(&fes);
    b.AddDomainIntegrator(new DomainLFIntegrator(zero));
    //b.Assemble();

    //1D case
    
    Array<int> dir_attr(mesh->bdr_attributes.Max());
    dir_attr = 1; //all attributes to be dirichlet's in 1d case only

    // 2D case
/*
    Array<int> dir_attr(mesh->bdr_attributes.Max());
    Array<int> neu_attr(mesh->bdr_attributes.Max());

    dir_attr = 0;
    neu_attr = 0;

    dir_attr[0] = 1;
    dir_attr[2] = 1;

    neu_attr[1] = 1;
    neu_attr[3] = 1;
*/
    //cout << "boundary attr: " << endl;
   // for(int i=0;i<bdr_attr.Size();i++)
    //    cout << bdr_attr[i] << endl;

    Array<int> anode_mkr(mesh->bdr_attributes.Max());
    Array<int> cathode_mkr(mesh->bdr_attributes.Max());

    lap_1d_bc(anode_mkr, cathode_mkr);
    //lap_2d_bc(anode_mkr,cathode_mkr);
    
    u.ProjectBdrCoefficient(VCoef,anode_mkr);
    u.ProjectBdrCoefficient(zero,cathode_mkr);

    cout << "initialized u:" << endl;
    cout << u << endl;   

    //for Dirichlet's conditions
    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(dir_attr, ess_tdof_list);
  
    for(int i=0;i<ess_tdof_list.Size();i++)
        cout << ess_tdof_list[i] << endl;
    //cout << *ess_tdof_list << endl;

    //for Neumanns conditions (for 2D)
   // b.AddBoundaryIntegrator(new BoundaryLFIntegrator(epCoef),neu_attr);
    b.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list,u,b,A,X,B);
    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, X, 1, 500, 1e-12, 0.0);
    a.RecoverFEMSolution(X,b,u);

    Vector u_nodes;
    u.GetNodalValues(u_nodes,1);

    cout << "u final:"<< endl;

    for(int i=0;i < u_nodes.Size();i++)
        cout << u_nodes[i] << endl;

    //GradientGridFunctionCoefficient E(&u); //will have to work on

    GridFunction E(&fes);

    u.GetDerivative(1,0,E); //works for 1d. (comp,der_comp,&der) comp = 1 for scalar func, der_comp = 0,1,2 for x,y,z
   
    //Vector E_nodes;
   // E.GetNodalValues(E_nodes,1);

    //cout << "E field" << endl;
    //for(int i=0;i<E_nodes.Size();i++)
      //  cout << E_nodes[i] << endl;

    cout << "vector dim " << endl;
    cout << E.VectorDim() << endl;

    Vector E_vec = E.GetTrueVector();
    cout << "E_vec:" << endl;
    for(int i; i < E_vec.Size();i++)
        cout << E_vec[i] << endl;

    return E_vec;
}

//returns the correct boundary markers for the cathode and anode of the 1D problem
void lap_1d_bc(Array<int> &anode, Array<int> &cathode){

    anode = 0;
    cathode = 0;

    anode[0] = 1;
    cathode[1] = 1;

}

//returns the correct boundary markers for the cathode and anode of the 1D problem 
//TODO: add Neumann conditions
void lap_2d_bc(Array<int> &anode, Array<int> &cathode){

    anode = 0;
    cathode = 0;

    anode[0] = 1;
    cathode[2] = 1;

}




//define list of parameters

