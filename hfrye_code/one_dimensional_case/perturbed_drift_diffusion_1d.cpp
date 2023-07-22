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
#include "params.hpp"
#include "fe_evolution.hpp"

#ifndef MFEM_USE_SLEPC
#error This examples requires that MFEM is build with MFEM_USE_SLEPC=YES
#endif

using namespace std;
using namespace mfem;

//Computes E field with given mesh
GridFunction electric_potential(int, double, double, FiniteElementSpace&, bool);

//Checks for Physics in limits of drift-diffusion equation
//BCs are kept dirichlets homogenous 
void time_indep_diffusion(FiniteElementSpace&, BlockMatrix&,Array<int> &);//elliptic problem
void time_dep_diffusion(FiniteElementSpace&, BlockMatrix&,BlockMatrix&, Array<int> &);//parabolic problem
void advection_dominated_flow(FiniteElementSpace&, BlockMatrix&, BlockMatrix&, Array<int> &);//

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

if(dim == 1){
   v(0) = 1.0;
} else if(dim == 2){
    v(0) = 1.0;
    v(1) = 1.0;
} else {
    cout << "dimension size not implmeneted yet. Aborting..."<<endl;
}
}

double time_dep_diff_ic(const Vector &x){
    int dim = x.Size();
    Vector X(dim);

    if(dim == 1){
        return 1.0;
    } else if(dim == 2){
        //v(0) = sin(M_PI*X(0));
        //v(1) = sin(M_PI*X(1));
        //return  sin(M_PI*x(0)) + sin(M_PI*x(1));
        //return 0.0;
        return exp((-pow(x(0)-0.5,2)/0.1)-(pow(x(1)-0.5,2)/0.1));
    } else {
    cout << "dimension size not implmeneted yet. Aborting..."<<endl;
    return 1;
    }
}
/*
double noise_ic(const Vector &x){
    int dim = x.Size();
    Vector X(dim);


}
*/
int main(int agrc, char *argv[]){
    const char *mesh_file = "../../data/inline-quad.mesh"; //Currently 2-d basic mesh to easily view results in visit
    //const char *mesh_file = "../../data/inline-segment.mesh";
    int order = 2; //second order legendre-Gauss solver (Make sure that's what it uses as according to Shibata paper)

    //Define Mesh 
    Mesh mesh(mesh_file);

    //Mesh mesh(10,1.0);
    // possibly refine mesh to increase resolution
    mesh.UniformRefinement();
   // mesh.UniformRefinement();

    int dim = mesh.Dimension();

    //define finite element space; for now, try to stick to H1 space and implement nonlinNSint but might need different space
    H1_FECollection fec(order,dim); //Gauss-Legendre unsupported
    FiniteElementSpace fespace(&mesh,&fec);
    cout << "Number of unknowns: " << fespace.GetTrueVSize() << endl;
    cout << "Vector Dim: " << fespace.GetVDim() << endl;

    //boundary markers for cathode and anode
    //For 1D case
    Array<int> cathode_bdr(mesh.bdr_attributes.Max());
    Array<int> anode_bdr(mesh.bdr_attributes.Max());
    cathode_bdr =0;
    anode_bdr = 0;
    cathode_bdr[0]=1;
    anode_bdr[2]=1; //change for 2d
    //E_field for 1d

    GridFunction e_pot = electric_potential(order,V,epsilon,fespace,true);
    GradientGridFunctionCoefficient E(&e_pot);
    
    ScalarVectorProductCoefficient v_e(mu_e,E);
    ScalarVectorProductCoefficient v_p(mu_p,E);

    FiniteElementSpace Vfespace(&mesh,&fec,dim);

    GridFunction v_e_grid(&Vfespace);
    v_e_grid.ProjectCoefficient(v_e);
    GridFunction v_p_grid(&Vfespace);
    v_p_grid.ProjectCoefficient(v_p);

    DivergenceGridFunctionCoefficient div_v_e(&v_e_grid);
    DivergenceGridFunctionCoefficient div_v_p(&v_p_grid);

    //Computing L2 norms of velocities
    InnerProductCoefficient v_e_mag_sqr(v_e,v_e);
    PowerCoefficient v_e_mag(v_e_mag_sqr,0.5);

   // InnerProductCoefficient v_p_mag_sqr(v_p,v_p);
    //PowerCoefficient v_p_mag(v_p_mag_sqr,0.5);

    SumCoefficient alpha_min_eta(Alpha,Eta,1.0,-1.0);
    ProductCoefficient ame_ve_mag(alpha_min_eta,v_e_mag);//alpha minus eta times magnitude v_e

    ProductCoefficient alpha_ve_mag(Alpha,v_e_mag);

    // constants to be multiplied into mass integrators for S matrices
    SumCoefficient SeeTot(ame_ve_mag,div_v_e,1.0,-1.0);
    SumCoefficient SpeTot(alpha_ve_mag,div_v_p,1.0,-1.0);

    //VectorFunctionCoefficient vel_test(dim,velocity_function);
    
    //BCs will have to specify in function the specifics later
    
    //For initial conditions primarily for physics check, will have to specify later
    
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

    M.AddDomainIntegrator(new MassIntegrator); //same for both equations
   
    Ae.AddDomainIntegrator(new AdvectionSUPGIntegrator(v_e,1.0,diff_const_e));//note that SUPG tau is dependant on diffusion
    De.AddDomainIntegrator(new DiffusionIntegrator(diff_const_e));   

    Kee.AddBoundaryIntegrator(new VectorNormedMassIntegrator(v_e,one), cathode_bdr); //TODO: will have to verify correct boundary term ;will need to split up into Kee and Kpe terms
    //Kee.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(velocity_e), cathode_bdr);

    Kep.AddBoundaryIntegrator(new VectorNormedMassIntegrator(v_p,Gamma), cathode_bdr);
    //Kep.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(velocity_p), cathode_bdr);
    See.AddDomainIntegrator(new MassIntegrator(SeeTot));
    Ap.AddDomainIntegrator(new AdvectionSUPGIntegrator(v_p,1.0,diff_const_p));
    Dp.AddDomainIntegrator(new DiffusionIntegrator(diff_const_p));
    Kp.AddBoundaryIntegrator(new VectorNormedMassIntegrator(v_p,one),anode_bdr);
    //Kep.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(velocity_p), anode_bdr);
    Spe.AddDomainIntegrator(new MassIntegrator(SpeTot));
    
    //after integrators are correctly defined, the Jacobian can be simply added or subtracted together 
    // as specified after assembling
    
    M.Assemble(); 
    Ae.Assemble(); 
    De.Assemble();   
    Kee.Assemble();
    Kep.Assemble();
    See.Assemble();

    Ap.Assemble();
    Dp.Assemble();
    Kp.Assemble();
    Spe.Assemble();
    // set up Jacobians

    SparseMatrix J11;

    //J11 = -Ae.SpMat() - De.SpMat() + Kee.SpMat() + See.SpMat();
    J11 = Ae.SpMat();
    J11 *= -1.0;
    J11.Add(-1.0,De.SpMat());
    J11.Add(1.0,Kee.SpMat());
    J11.Add(1.0,See.SpMat());

    SparseMatrix J12; 
    J12 = Kep.SpMat();

 cout <<  "J12" << J12.GetData() << endl;

    SparseMatrix J21; 
    J21 = Spe.SpMat();

 cout <<  "J21" << J21.GetData() << endl;
    SparseMatrix J22;

    //J22 = -Ap.SpMat() - Dp.SpMat() + Kp.SpMat();
    J22 = Ap.SpMat();
    J22 *= -1.0;
    J22.Add(-1.0,Dp.SpMat());
    J22.Add(1.0,Kp.SpMat());

 cout <<  "J22" << J22.GetData() << endl;

    //jacobian block check
    cout<< "Size check: J11 "<<J11.Size()<< " J12 "<<J12.Size()<< " J21 "<<J21.Size()<< " J22 "<<J22.Size()<<endl;

    //Set up block jacobian and block mass matrix
    Array<int> block_offsets(3); // number of variables + 1
    for(int k = 0; k < 3; k++){
        block_offsets[k] = k * fespace.GetNDofs();
    }

    BlockMatrix Block_Jacobian(block_offsets);
    Block_Jacobian.SetBlock(0,0,&J11);
    Block_Jacobian.SetBlock(0,1,&J12);
    Block_Jacobian.SetBlock(1,0,&J21);
    Block_Jacobian.SetBlock(1,1,&J22);

    Block_Jacobian.Finalize();

   //for(int i =0;i<J11.Size()+J21.Size();i++)
    //    cout << "Block Jacobian Diag: "<<Block_Jacobian.Elem(i,i)<<endl;

    SparseMatrix MassMat(M.SpMat());

    BlockMatrix Block_Mass(block_offsets);
    Block_Mass.SetBlock(0,0,&MassMat);
    Block_Mass.SetBlock(1,1,&MassMat);

    Block_Mass.Finalize();

    //from here the eigenvalue problem can be solved

    //physics tests

    time_indep_diffusion(fespace, Block_Jacobian, block_offsets);//elliptic problem
    time_dep_diffusion(fespace, Block_Jacobian, Block_Mass, block_offsets);//parabolic problem

    //eigensolver slepc
    int nev = 5;
    SlepcEigenSolver *slepc = new SlepcEigenSolver();
    slepc->SetNumModes(nev);
    slepc->SetWhichEigenpairs(SlepcEigenSolver::Target_REAL);
    slepc->SetTarget(0.0);
    slepc->SetSpectralTransformation(SlepcEigenSolver::SHIFT_INVERT);
    slepc->SetOperators(Block_Jacobian,Block_Mass); //might need to convert datatypes to petsc matrices

    Array<double> eig_vals;

    slepc->Solve();
    eig_vals.SetSize(nev);
    for(int i = 0; i<nev; i++){
        slepc->GetEigenvalue(i,eig_vals[i]);
    }

    BlockVector u_block(block_offsets); //to intialize size of matrix containing eigenvectors
    Matrix eig_vecs(u_block.Size(),nev); //eigenvectors stored as columns in this matrix

    for(int i = 0; i < nev; i++){
        slepc->GetEigenvector(i,u_block);
        for(int j = 0; j < u_block.Size(); j++){
            eig_vecs(j,i) = u_block[j];
        }
    }
    
/*
    DataCollection *dc = NULL;
    dc = new VisItDataCollection("Eigenmodes_test",&mesh);
    dc->SetPrecision(8);
    for()
    dc->RegisterField("solution",&u);
    dc->Save();
    delete dc;
*/
   //ofstream meshfile;
    //meshfile.open("mesh.dat");
    //for(int i=0; i < mesh.GetElementSize)

    return 0;
}

//solves preliminary laplace problem for the electric potential and then returns the resulting electric field
GridFunction electric_potential(int order, double V, double epsilon, FiniteElementSpace &fes, bool print_epot){// pass in boundary elements along with fes since it should match up with d-d problem 
//permittivity
    //Mesh &mesh = fes.GetMesh();
    Mesh &mesh = *fes.GetMesh();
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

    Array<int> anode_mkr(mesh.bdr_attributes.Max());
    Array<int> cathode_mkr(mesh.bdr_attributes.Max());

    cout << "mesh dimension: " << mesh.Dimension() << endl;
    int dim = 2;

    //if(dim==1){
    //1D case
/*
    Array<int> dir_attr(mesh.bdr_attributes.Max());
    dir_attr = 1; //all attributes to be dirichlet's in 1d case only

    lap_1d_bc(anode_mkr, cathode_mkr);
*/
   // } else if(dim==2){
    // 2D case

    Array<int> dir_attr(mesh.bdr_attributes.Max());
    //Array<int> neu_attr(mesh.bdr_attributes.Max());

    dir_attr = 0;
    //neu_attr = 0;

    dir_attr[0] = 1;
    dir_attr[2] = 1;

    lap_2d_bc(anode_mkr,cathode_mkr);
  //  }
/*
    neu_attr[1] = 1;
    neu_attr[3] = 1;
*/
    //cout << "boundary attr: " << endl;
   // for(int i=0;i<bdr_attr.Size();i++)
    //    cout << bdr_attr[i] << endl;
    
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

    cout << "u final nodal:"<< endl;

    for(int i=0;i < u_nodes.Size();i++)
        cout << u_nodes[i] << endl;

     // for diagnostic purposes
     if(print_epot)
     {
        DataCollection *dc = NULL;
        dc = new VisItDataCollection("Electric_Potential",&mesh);
        dc->SetPrecision(8);
        dc->RegisterField("solution",&u);
        dc->Save();
        delete dc;
     }

     cout << "GetNE check : " << mesh.GetNE() << endl;

    return u;
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

//Where time derivative in equations go to zero, and velocity is set to zero
void time_indep_diffusion(FiniteElementSpace &fes, BlockMatrix& J, Array<int> &offsets){
    Mesh &mesh = *fes.GetMesh();
    GridFunction n_e(&fes), n_p(&fes);

    //MemoryType mt = device.GetMemoryType();
    BlockVector n_block(offsets), rhs(offsets);
    //n_block = 0;

    ConstantCoefficient one(1.0);
    LinearForm e_source(&fes);
    e_source.AddDomainIntegrator(new DomainLFIntegrator(one));
    e_source.Assemble();

    LinearForm p_source(&fes);
    p_source.AddDomainIntegrator(new DomainLFIntegrator(one));
    p_source.Assemble();

    //init rhs of test problem
    rhs = 0.0;
    rhs.GetBlock(0) = e_source;
    rhs.GetBlock(1) = p_source;

    //init solutions
    n_block.GetBlock(0) = n_e;
    n_block.GetBlock(1) = n_p;

    //ensuring homogenous Dirichlet boundary conditions
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1;

    SparseMatrix *JSp = J.CreateMonolithic();
    GSSmoother M(*JSp);
    GMRES(*JSp,M,rhs,n_block, 1, 500, 10, 1e-12, 0.0);

    n_e.MakeRef(&fes,n_block.GetBlock(0),0);
    n_p.MakeRef(&fes,n_block.GetBlock(1),0);

    DataCollection *dc = NULL;
    dc = new VisItDataCollection("Time_ind_diff",&mesh);
    dc->SetPrecision(8);
    dc->RegisterField("solution_e",&n_e);
    dc->RegisterField("solution_p",&n_p);
    dc->Save();
    delete dc;

}

void time_dep_diffusion(FiniteElementSpace& fes, BlockMatrix& J, BlockMatrix& M, Array<int> &offsets){
    Mesh &mesh = *fes.GetMesh();
    GridFunction n_e(&fes), n_p(&fes);

    int dim = mesh.Dimension();

    SparseMatrix *JSp = J.CreateMonolithic();
    SparseMatrix *MSp = M.CreateMonolithic();

    ODESolver *ode_solver = new BackwardEulerSolver;

    FunctionCoefficient ne_0(time_dep_diff_ic);
    FunctionCoefficient np_0(time_dep_diff_ic);

    n_e.ProjectCoefficient(ne_0);
    n_p.ProjectCoefficient(np_0);

    //TODO:define block vector
    BlockVector n_block(offsets), b(offsets);
    n_block.GetBlock(0) = n_e;
    n_block.GetBlock(1) = n_p;

    LinearForm b_e(&fes);
    b_e.AddDomainIntegrator(new DomainLFIntegrator(zero));
    b_e.Assemble();

    LinearForm b_p(&fes);
    b_p.AddDomainIntegrator(new DomainLFIntegrator(zero));
    b_p.Assemble();

    b.GetBlock(0) = b_e;
    b.GetBlock(1) = b_p;
//initial conditions
    DataCollection *dc = NULL;
    dc = new VisItDataCollection("Time_dep_problem_1V",&mesh);
    dc->SetPrecision(8);
    dc->RegisterField("solution_e",&n_e);
    dc->RegisterField("solution_p",&n_p);
    dc->Save();

    cout << "n_e init:"<<endl;
    cout<< n_e << endl;

    FE_Evolution diff(*MSp,*JSp,b);

    double t = 0.0;
    diff.SetTime(t);
    ode_solver->Init(diff);

    int vis_steps = 5;
    bool done = false;
    for (int ti = 0; !done; )
    {
        
        double dt_real = min(dt, t_final - t);
        ode_solver->Step(n_block, t, dt_real);
        ti++;

        n_e.MakeRef(&fes,n_block.GetBlock(0),0);
        n_p.MakeRef(&fes,n_block.GetBlock(1),0);

        //cout << "n_e check" << n_e << endl;

        done = (t >= t_final - 1e-8*dt);

        if (done || ti % vis_steps == 0)
        {
            cout << "time step: " << ti << ", time: " << t << endl;

                dc->SetCycle(ti);
                dc->SetTime(t);
                dc->Save();
        }
        
    }

    n_e.MakeRef(&fes,n_block.GetBlock(0),0);
    n_p.MakeRef(&fes,n_block.GetBlock(1),0);

    cout << "n_e final:"<<endl;
    cout<< n_e << endl;

    dc = NULL;
    dc = new VisItDataCollection("Time_dep_final_1V",&mesh);
    dc->SetPrecision(8);
    dc->RegisterField("solution_e",&n_e);
    dc->RegisterField("solution_p",&n_p);
    dc->Save();
    
    delete ode_solver;
    delete dc;

}


