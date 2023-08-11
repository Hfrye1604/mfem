/*1-d plate to plate conduction discharge tube code using MFEM for the Drift-Diffusion equations
    --The purpose of this code is to attempt to solve the D-D perturbation equations
    to verify the physics before defining for the jacobian and mass matrices of the system
    for stability analysis. This code is written with dimensional generality in mind.
*/
// #include "/g/g11/frye11/mfem/mfem.hpp"
#include "../../mfem.hpp"
// #include "NSnonlininteg_modified.hpp"
// #include "fem/advectionSUPGinteg.hpp"
#include <fstream>
#include <iostream>
#include <math.h>
#include "params.hpp"
#include "fe_evolution.hpp"

#ifndef MFEM_USE_SLEPC
#error This examples requires that MFEM is build with MFEM_USE_SLEPC=YES
#endif

using namespace std;
using namespace mfem;

// Computes E field with given mesh
GridFunction electric_potential(int, double, double, FiniteElementSpace &, bool);

// Checks for Physics in limits of drift-diffusion equation
// BCs are kept dirichlets homogenous
void time_indep_diffusion(FiniteElementSpace &, BlockMatrix &, Array<int> &, bool);              // elliptic problem
void time_dep_diffusion(FiniteElementSpace &, BlockMatrix &, BlockMatrix &, Array<int> &, bool); // parabolic problem
// void advection_dominated_flow(FiniteElementSpace&, BlockMatrix&, BlockMatrix&, Array<int> &);

int main(int argc, char *argv[])
{

    // Define Mesh
    Mesh mesh(mesh_file);

    // Mesh mesh(10,1.0);
    //  possibly refine mesh to increase resolution
     mesh.UniformRefinement();
     //mesh.UniformRefinement();
     //mesh.UniformRefinement();
    //mesh.UniformRefinement();

    int dim = mesh.Dimension();

    // define finite element space; for now, try to stick to H1 space and implement nonlinNSint but might need different space
    H1_FECollection fec(order, dim); // Gauss-Legendre unsupported
    FiniteElementSpace fespace(&mesh, &fec);
    cout << "Number of unknowns: " << fespace.GetTrueVSize() << endl;
    cout << "Vector Dim: " << fespace.GetVDim() << endl;

    // boundary markers for cathode and anode
    // For 2D parallel plates
    Array<int> bdr_mkr(mesh.bdr_attributes.Max()); // for essential BCs
    Array<int> cathode_bdr(mesh.bdr_attributes.Max());
    Array<int> anode_bdr(mesh.bdr_attributes.Max());

     lap_2d_bc(anode_bdr,cathode_bdr, bdr_mkr);
    //lap_sqr_disc(anode_bdr, cathode_bdr, bdr_mkr);
    // l_shape_bc(anode_bdr, cathode_bdr, bdr_mkr);
     //corner_plane_bc(anode_bdr, cathode_bdr, bdr_mkr);
    //lap_1d_bc(anode_bdr, cathode_bdr, bdr_mkr);

    GridFunction e_pot = electric_potential(order, V, epsilon, fespace, true);
    GradientGridFunctionCoefficient E(&e_pot);

    InnerProductCoefficient E_mag_sqr(E, E);
    PowerCoefficient E_mag(E_mag_sqr, 0.5);

    // mu_e
    TransformedCoefficient mu_e(&E_mag, piecewise_func);

    
    // diff_const_e
    ConstantCoefficient one_over_p_100(1.0 / (pressure * 100));
    ProductCoefficient E_over_100p(E_mag, one_over_p_100);
    PowerCoefficient E_o_100p_powered_2807(E_over_100p, .2807);
    SumCoefficient sum1(E_o_100p_powered_2807, E_over_100p, 22.26, 0.5787);
    SumCoefficient diff_const_e(46.0, sum1, 1.0, 1.0);

   
    // mu_p
    SumCoefficient exp_term(1251, E_over_100p, -0.005, -0.005);
    TransformedCoefficient exp_E_term(&exp_term, exp);
    SumCoefficient mu_p(.000046, exp_E_term, 1.0, 0.0801);

    
    // diff_const_p
    PowerCoefficient E_o_100p_powered_8378(E_over_100p, 0.8378);
    SumCoefficient diff_const_p(0.645163, E_o_100p_powered_8378, 1.0, 0.0562447);

    // ProductCoefficient Alpha = Alpha_func(E);
    // Alpha
    double A_param = 2900;
    double B_param = 28070; // expressions uses from table 14.1 in Lieberman paper found through schynder paper
    PowerCoefficient one_over_E_mag(E_mag_sqr, -0.5);
    ConstantCoefficient nB_p(-B_param * pressure);
    ProductCoefficient nB_p_over_E(nB_p, one_over_E_mag);
    TransformedCoefficient exponential_coeff(&nB_p_over_E, exp);
    ConstantCoefficient p_A(pressure * A_param);
    ProductCoefficient Alpha(p_A, exponential_coeff);

    // TransformedCoefficient Eta = Eta_func(E);
    // Eta
    TransformedCoefficient Eta(&E_mag, interpolate_eta);

    // velocity vector fields
    ScalarVectorProductCoefficient v_e(mu_e, E);
    ScalarVectorProductCoefficient v_p(mu_p, E);

    FiniteElementSpace Vfespace(&mesh, &fec, dim);

    cout << "v_e check: " << v_e.GetVDim() << endl;
    cout << "dim check: " << dim << endl;

    // velocity field divergences
    GridFunction v_e_grid(&Vfespace);
    v_e_grid.ProjectCoefficient(v_e);
    GridFunction v_p_grid(&Vfespace);
    v_p_grid.ProjectCoefficient(v_p);

    DivergenceGridFunctionCoefficient div_v_e(&v_e_grid);
    DivergenceGridFunctionCoefficient div_v_p(&v_p_grid);

    // Computing L2 norms of velocities
    InnerProductCoefficient v_e_mag_sqr(v_e, v_e);
    PowerCoefficient v_e_mag(v_e_mag_sqr, 0.5);

    // InnerProductCoefficient v_p_mag_sqr(v_p,v_p);
    // PowerCoefficient v_p_mag(v_p_mag_sqr,0.5);

    SumCoefficient alpha_min_eta(Alpha, Eta, 1.0, -1.0);
    ProductCoefficient ame_ve_mag(alpha_min_eta, v_e_mag); // alpha minus eta times magnitude v_e

    ProductCoefficient alpha_ve_mag(Alpha, v_e_mag);

    // constants to be multiplied into mass integrators for S matrices
    SumCoefficient SeeTot(ame_ve_mag, div_v_e, 1.0, -1.0);
    SumCoefficient SpeTot(alpha_ve_mag, div_v_p, 1.0, -1.0);

    // coefficients check
    GridFunction mu_e_grid(&fespace);
    mu_e_grid.ProjectCoefficient(mu_e);
    GridFunction diff_const_e_grid(&fespace);
    diff_const_e_grid.ProjectCoefficient(diff_const_e);
    GridFunction mu_p_grid(&fespace);
    mu_p_grid.ProjectCoefficient(mu_p);
    GridFunction diff_const_p_grid(&fespace);
    diff_const_p_grid.ProjectCoefficient(diff_const_p);

    GridFunction Alpha_grid(&fespace);
    Alpha_grid.ProjectCoefficient(Alpha);
    GridFunction Eta_grid(&fespace);
    Eta_grid.ProjectCoefficient(Eta);

    //E_mag
    GridFunction E_mag_grid(&fespace);
    E_mag_grid.ProjectCoefficient(E_mag);

    DataCollection *dc1 = NULL;
    dc1 = new VisItDataCollection("Coefficients", &mesh);
    dc1->SetPrecision(8);
    dc1->RegisterField("mu_e", &mu_e_grid);
    dc1->RegisterField("diff_const_e", &diff_const_e_grid);
    dc1->RegisterField("mu_p", &mu_p_grid);
    dc1->RegisterField("diff_const_p", &diff_const_p_grid);
    dc1->RegisterField("Alpha", &Alpha_grid);
    dc1->RegisterField("Eta", &Eta_grid);
    dc1->RegisterField("|E|",&E_mag_grid);
    dc1->Save();
    delete dc1;

    // BCs will have to specify in function the specifics later

    // For initial conditions primarily for physics check, will have to specify later

    // Following integrator matrix notation from Shibata paper
    BilinearForm M(&fespace); // mass matrix the same for both species
    BilinearForm Ae(&fespace);
    BilinearForm De(&fespace);
    BilinearForm Kee(&fespace);
    BilinearForm Kep(&fespace);
    BilinearForm See(&fespace);

    BilinearForm Ap(&fespace);
    BilinearForm Dp(&fespace);
    BilinearForm Kp(&fespace);
    BilinearForm Spe(&fespace);

    M.AddDomainIntegrator(new MassIntegrator); // same for both equations
/*
    Ae.AddDomainIntegrator(new AdvectionSUPGIntegrator(v_e, 1.0, diff_const_e)); // note that SUPG tau is dependant on diffusion
    De.AddDomainIntegrator(new DiffusionIntegrator(one));
    Kee.AddBoundaryIntegrator(new VectorNormedMassIntegrator(v_e, one), cathode_bdr); // TODO: will have to verify correct boundary term ;will need to split up into Kee and Kpe terms
    Kep.AddBoundaryIntegrator(new VectorNormedMassIntegrator(v_p, Gamma), cathode_bdr);
    See.AddDomainIntegrator(new MassIntegrator(one));

    Ap.AddDomainIntegrator(new AdvectionSUPGIntegrator(v_p, 1.0, diff_const_p));
    Dp.AddDomainIntegrator(new DiffusionIntegrator(one));
    Kp.AddBoundaryIntegrator(new VectorNormedMassIntegrator(v_p, one), anode_bdr);
    Spe.AddDomainIntegrator(new MassIntegrator(one));
*/

    Ae.AddDomainIntegrator(new AdvectionSUPGIntegrator(v_e, 1.0, diff_const_e)); // note that SUPG tau is dependant on diffusion
    De.AddDomainIntegrator(new DiffusionIntegrator(diff_const_e));
    Kee.AddBoundaryIntegrator(new VectorNormedMassIntegrator(v_e, one), cathode_bdr); // TODO: will have to verify correct boundary term ;will need to split up into Kee and Kpe terms
    Kep.AddBoundaryIntegrator(new VectorNormedMassIntegrator(v_p, Gamma), cathode_bdr);
    See.AddDomainIntegrator(new MassIntegrator(SeeTot));

    Ap.AddDomainIntegrator(new AdvectionSUPGIntegrator(v_p, 1.0, diff_const_p));
    Dp.AddDomainIntegrator(new DiffusionIntegrator(diff_const_p));
    Kp.AddBoundaryIntegrator(new VectorNormedMassIntegrator(v_p, one), anode_bdr);
    Spe.AddDomainIntegrator(new MassIntegrator(SpeTot));


    // after integrators are correctly defined, the Jacobian can be simply added or subtracted together
    //  as specified after assembling

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

    // J11 = -Ae.SpMat() - De.SpMat() + Kee.SpMat() + See.SpMat();
    J11 = Ae.SpMat();
    J11 *= -1.0;
    J11.Add(-1.0, De.SpMat());
    J11.Add(1.0, Kee.SpMat());
    J11.Add(1.0, See.SpMat());

    cout << "J11" << J11.GetData() << endl;

    SparseMatrix J12;
    J12 = Kep.SpMat();

    cout << "J12" << J12.GetData() << endl;

    SparseMatrix J21;
    J21 = Spe.SpMat();

    cout << "J21" << J21.GetData() << endl;
    SparseMatrix J22;

    // J22 = -Ap.SpMat() - Dp.SpMat() + Kp.SpMat();
    J22 = Ap.SpMat();
    J22 *= -1.0;
    J22.Add(-1.0, Dp.SpMat());
    J22.Add(1.0, Kp.SpMat());

    cout << "J22" << J22.GetData() << endl;

    // jacobian block check
    cout << "Size check: J11 " << J11.Size() << " J12 " << J12.Size() << " J21 " << J21.Size() << " J22 " << J22.Size() << endl;

    // Set up block jacobian and block mass matrix
    Array<int> block_offsets(3); // number of variables + 1
    for (int k = 0; k < 3; k++)
    {
        block_offsets[k] = k * fespace.GetNDofs();
    }

    BlockMatrix Block_Jacobian(block_offsets);
    Block_Jacobian.SetBlock(0, 0, &J11);
    Block_Jacobian.SetBlock(0, 1, &J12);
    Block_Jacobian.SetBlock(1, 0, &J21);
    Block_Jacobian.SetBlock(1, 1, &J22);

    Block_Jacobian.Finalize();

    // for(int i =0;i<J11.Size()+J21.Size();i++)
    //    cout << "Block Jacobian Diag: "<<Block_Jacobian.Elem(i,i)<<endl;
    // for(int i = 0; i < J11.Size();i++){
    // for(int j = 0; j < J11.Size();j++){
    //   if(J11(j,i)){
    //      cout << Block_Jacobian(j,i) << " ";
    //  }
    //}
    // cout << endl;
    //}

    SparseMatrix MassMat(M.SpMat());

    BlockMatrix Block_Mass(block_offsets);
    Block_Mass.SetBlock(0, 0, &MassMat);
    Block_Mass.SetBlock(1, 1, &MassMat);

    Block_Mass.Finalize();

    // from here the eigenvalue problem can be solved

    // physics tests
   
    time_indep_diffusion(fespace, Block_Jacobian, block_offsets, run_time_indep);             // elliptic problem
    time_dep_diffusion(fespace, Block_Jacobian, Block_Mass, block_offsets, run_fe_evolution); // parabolic problem
   
    // eigensolver slepc

    Mpi::Init(argc, argv);
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

    // first need to convert from sparseMatrix to PetscParMatrix
    MFEMInitializeSlepc(NULL, NULL, NULL, NULL);
    //PetscComplex imag = PETSC_I
    //SparseMatrix test();

    PetscParMatrix *pJ = new PetscParMatrix(Block_Jacobian.CreateMonolithic());
    PetscParMatrix *pM = new PetscParMatrix(Block_Mass.CreateMonolithic());

    int nev = 10;
    SlepcEigenSolver *slepc = new SlepcEigenSolver(MPI_COMM_WORLD);
    slepc->SetNumModes(nev);
    slepc->SetWhichEigenpairs(SlepcEigenSolver::TARGET_REAL);
    slepc->SetTarget(0.0);
    slepc->SetSpectralTransformation(SlepcEigenSolver::SHIFT_INVERT);
    slepc->SetOperators(*pJ, *pM); // might need to convert datatypes to petsc

    // slepc->SetTol(1e-6);
    // slepc->SetMaxIter(1500);

    Array<double> eig_vals;

    slepc->Solve();

    cout << "num converged: " << slepc->GetNumConverged() << endl;

    eig_vals.SetSize(nev);

    for (int i = 0; i < nev; i++)
    {
        slepc->GetEigenvalue(i, eig_vals[i]);
    }
    
    BlockVector u_block(block_offsets);        // to intialize size of matrix containing eigenvectors
    DenseMatrix eig_vecs(u_block.Size(), nev); // eigenvectors stored as columns in this matrix

    for (int i = 0; i < nev; i++)
    {
        slepc->GetEigenvector(i, u_block);
        for (int j = 0; j < u_block.Size(); j++)
        {
            eig_vecs(j, i) = u_block[j];
        }
    }

    BlockVector first_eig_vec(eig_vecs.GetColumn(0), block_offsets);
    GridFunction ne_first_eigenmode(&fespace);
    GridFunction np_first_eigenmode(&fespace);

    BlockVector second_eig_vec(eig_vecs.GetColumn(1), block_offsets);
    GridFunction ne_second_eigenmode(&fespace);
    GridFunction np_second_eigenmode(&fespace);

    BlockVector third_eig_vec(eig_vecs.GetColumn(2), block_offsets);
    GridFunction ne_third_eigenmode(&fespace);
    GridFunction np_third_eigenmode(&fespace);

    ne_first_eigenmode.MakeRef(&fespace, first_eig_vec.GetBlock(0), 0);
    ne_second_eigenmode.MakeRef(&fespace, second_eig_vec.GetBlock(0), 0);
    ne_third_eigenmode.MakeRef(&fespace, third_eig_vec.GetBlock(0), 0);

    np_first_eigenmode.MakeRef(&fespace, first_eig_vec.GetBlock(1), 0);
    np_second_eigenmode.MakeRef(&fespace, second_eig_vec.GetBlock(1), 0);
    np_third_eigenmode.MakeRef(&fespace, third_eig_vec.GetBlock(1), 0);

    cout << "Eigenvalues:" << endl;
    for (int i = 0; i < nev; i++)
        cout << eig_vals[i] << endl;

    DataCollection *dc = NULL;
    dc = new VisItDataCollection("Eigenmodes", &mesh);
    dc->SetPrecision(8);
    dc->RegisterField("n_e 1st_eigenmode", &ne_first_eigenmode);
    dc->RegisterField("n_e 2nd_eigenmode", &ne_second_eigenmode);
    dc->RegisterField("n_e 3rd_eigenmode", &ne_third_eigenmode);
    dc->RegisterField("n_p 1st_eigenmode", &np_first_eigenmode);
    dc->RegisterField("n_p 2nd_eigenmode", &np_second_eigenmode);
    dc->RegisterField("n_p 3rd_eigenmode", &np_third_eigenmode);
    dc->Save();
    delete dc;


    if(dim == 1){ //for use in matlab
        GridFunction segment(&fespace);
        Vector seg_vec;
        Vector n_e_eigmode1;
        Vector n_p_eigmode1;
        mesh.GetNodes(segment);
        segment.GetNodalValues(seg_vec,1);
        ne_first_eigenmode.GetNodalValues(n_e_eigmode1,1);
        np_first_eigenmode.GetNodalValues(n_p_eigmode1,1);

        for(int i=0; i < seg_vec.Size(); i++)
            cout << seg_vec.Elem(i) << " ";
        cout<<endl;\
        cout<<endl;
        for(int i=0; i < seg_vec.Size(); i++)
            cout << n_e_eigmode1.Elem(i) << " ";
        cout<<endl;
        cout<<endl;
        for(int i=0; i < seg_vec.Size(); i++)
            cout << n_p_eigmode1.Elem(i) << " ";
   // dc->RegisterField("line_segment", mesh.GetNodes());
    }


    // ofstream meshfile;
    // meshfile.open("mesh.dat");
    // for(int i=0; i < mesh.GetElementSize)
    delete pJ;
    delete pM;
    delete slepc;

    return 0;
}

// solves preliminary laplace problem for the electric potential and then returns the resulting electric field
GridFunction electric_potential(int order, double V, double epsilon, FiniteElementSpace &fes, bool print_epot)
{   // pass in boundary elements along with fes since it should match up with d-d problem
    // permittivity
    // Mesh &mesh = fes.GetMesh();
    Mesh &mesh = *fes.GetMesh();
    // Array<int> ess_bdr(mesh.bdr_attributes.Max());
    //  ess_bdr = 1;
    // fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // u represents the electric potential
    GridFunction u(&fes);
    u = 0.0;
    // u[0]=V;

    cout.precision(8);
    cout << "initialized :" << endl;
    cout << u << endl;

    ConstantCoefficient VCoef(V); // anode voltage
    ConstantCoefficient epCoef(epsilon);

    BilinearForm a(&fes);
    a.AddDomainIntegrator(new DiffusionIntegrator(epCoef));
    a.Assemble();

    ConstantCoefficient zero(0.0);
    LinearForm b(&fes);
    b.AddDomainIntegrator(new DomainLFIntegrator(zero));
    // b.Assemble();

    Array<int> anode_bdr(mesh.bdr_attributes.Max());
    Array<int> cathode_bdr(mesh.bdr_attributes.Max());
    Array<int> bdr_mkr(mesh.bdr_attributes.Max());

    cout << "mesh dimension: " << mesh.Dimension() << endl;

   // lap_1d_bc(anode_bdr,cathode_bdr,bdr_mkr);
     lap_2d_bc(anode_bdr,cathode_bdr,bdr_mkr);
    // lap_2d_l_shape(anode_mkr, cathode_mkr);
    //lap_sqr_disc(anode_bdr, cathode_bdr, bdr_mkr);
    // l_shape_bc(anode_mkr,cathode_mkr,dir_attr);
     //corner_plane_bc(anode_bdr, cathode_bdr, bdr_mkr);

    //  }
    /*
        neu_attr[1] = 1;
        neu_attr[3] = 1;
    */
    // cout << "boundary attr: " << endl;
    // for(int i=0;i<bdr_attr.Size();i++)
    //    cout << bdr_attr[i] << endl;

    u.ProjectBdrCoefficient(VCoef, anode_bdr);
    u.ProjectBdrCoefficient(zero, cathode_bdr);

    cout << "initialized u:" << endl;
    cout << u << endl;

    // for Dirichlet's conditions
    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(bdr_mkr, ess_tdof_list);

    for (int i = 0; i < ess_tdof_list.Size(); i++)
        cout << ess_tdof_list[i] << endl;
    // cout << *ess_tdof_list << endl;

    // for Neumanns conditions (for 2D)
    // b.AddBoundaryIntegrator(new BoundaryLFIntegrator(epCoef),neu_attr);
    b.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, u, b, A, X, B);
    GSSmoother M((SparseMatrix &)(*A));
    PCG(*A, M, B, X, 1, 500, 1e-12, 0.0);
    a.RecoverFEMSolution(X, b, u);

    Vector u_nodes;
    u.GetNodalValues(u_nodes, 1);

    // cout << "u final nodal:"<< endl;

    // for(int i=0;i < u_nodes.Size();i++)
    //    cout << u_nodes[i] << endl;

    // for diagnostic purposes
    if (print_epot)
    {
        DataCollection *dc = NULL;
        dc = new VisItDataCollection("Electric_Potential", &mesh);
        dc->SetPrecision(8);
        dc->RegisterField("solution", &u);
        dc->Save();
        delete dc;
    }

    cout << "GetNE check : " << mesh.GetNE() << endl;

    return u;
}

// Where time derivative in equations go to zero, and velocity is set to zero
void time_indep_diffusion(FiniteElementSpace &fes, BlockMatrix &J, Array<int> &offsets, bool run)
{
    if (run == false)
        return;

    Mesh &mesh = *fes.GetMesh();
    GridFunction n_e(&fes), n_p(&fes);

    // MemoryType mt = device.GetMemoryType();
    BlockVector n_block(offsets), rhs(offsets);
    // n_block = 0;

    ConstantCoefficient one(1.0);
    LinearForm e_source(&fes);
    e_source.AddDomainIntegrator(new DomainLFIntegrator(one));
    e_source.Assemble();

    LinearForm p_source(&fes);
    p_source.AddDomainIntegrator(new DomainLFIntegrator(one));
    p_source.Assemble();

    // init rhs of test problem
    rhs = 0.0;
    rhs.GetBlock(0) = e_source;
    rhs.GetBlock(1) = p_source;

    // init solutions
    n_block.GetBlock(0) = n_e;
    n_block.GetBlock(1) = n_p;

    // ensuring homogenous Dirichlet boundary conditions
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1;

    SparseMatrix *JSp = J.CreateMonolithic();

    // iterative method
    // GSSmoother Prec(*JSp);
    // GMRES(*JSp,Prec,rhs,n_block, 1, 500, 10, 1e-12, 0.0);

    // direct
    DenseMatrix *JD = JSp->ToDenseMatrix();
    MatrixInverse *JDinv = JD->Inverse();
    JDinv->Mult(rhs, n_block);

    /*
    MFEMInitializePetsc();

    _p_Mat *pJ = new PetscParMatrix(J.CreateMonolithic());
    _p_Vec *pn = new PetscParVector(n_block);
    _p_Vec *prhs = new PetscParVector(rhs);

    _p_KSP *ksp;
    //_p_Vec *vec;
    KSPCreate(MPI_COMM_WORLD,ksp);
    KSPSetType(ksp,KSPNONE);
    KSPSetOperators(ksp,pJ);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp,prhs,pn);

    MFEMFinalizePetsc();
*/
    n_e.MakeRef(&fes, n_block.GetBlock(0), 0);
    n_p.MakeRef(&fes, n_block.GetBlock(1), 0);

    DataCollection *dc = NULL;
    dc = new VisItDataCollection("Time_ind_diff", &mesh);
    dc->SetPrecision(8);
    dc->RegisterField("solution_e", &n_e);
    dc->RegisterField("solution_p", &n_p);
    dc->Save();
    delete dc;

    // delete pJ;
}

void time_dep_diffusion(FiniteElementSpace &fes, BlockMatrix &J, BlockMatrix &M, Array<int> &offsets, bool run)
{

    if (run == false)
        return;
    cout << "time_dep_diff: Here?" << endl;
    Mesh &mesh = *fes.GetMesh();
    GridFunction n_e(&fes), n_p(&fes);

    int dim = mesh.Dimension();

    SparseMatrix *JSp_p = J.CreateMonolithic();
    SparseMatrix *MSp_p = M.CreateMonolithic();
    SparseMatrix JSp(*JSp_p);
    SparseMatrix MSp(*MSp_p);
    cout << "time_dep_diff: Here?" << endl;

    JSp.Finalize();
    MSp.Finalize();

    //  for(int i=0;i<JSp.Size();i++){
    // for(int j=0;j < JSp->Size();j++){
    //         cout << JSp.Elem(i,i) << endl;// " ";
    //}
    // cout << endl;
    // }

    ODESolver *ode_solver = new BackwardEulerSolver;
    cout << "time_dep_diff: Here?" << endl;

    // FunctionCoefficient ne_0(time_dep_diff_ic); Gaussian
    // FunctionCoefficient np_0(time_dep_diff_ic);

    FunctionCoefficient ne_0(noise_ic);
    FunctionCoefficient np_0(noise_ic);

    n_e.ProjectCoefficient(ne_0);
    n_p.ProjectCoefficient(np_0);

    // cout << "ne_0 check:" << endl;
    // cout << n_e << endl;

    // TODO:define block vector
    BlockVector n_block(offsets), b(offsets);
    // n_block = new BlockVector(offsets);
    // b = new BlockVector(offsets);
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
    // initial conditions
    DataCollection *dc = NULL;
    dc = new VisItDataCollection(FE_db_name, &mesh);
    dc->SetPrecision(8);
    dc->RegisterField("solution_e", &n_e);
    dc->RegisterField("solution_p", &n_p);
    dc->Save();

    // block preconditioner
    // cout << "n_e init:"<<endl;
    // cout<< n_e << endl;

    double t = 0.0;

    SparseMatrix A;
    Vector n_temp(n_block.Size());
    Vector n_next(n_block.Size());
    DenseMatrix *AD = NULL;

    if (!direct_solve)
    {
        FE_Evolution diff(MSp, JSp, b, offsets);
        diff.SetTime(t);
        ode_solver->Init(diff);
    }
    else
    {
        A = JSp;
        A *= -dt;
        A += MSp;
        AD = A.ToDenseMatrix();
        AD->Invert();
    }

    int vis_steps = 5;
    bool done = false;
    for (int ti = 0; !done;)
    {

        double dt_real = min(dt, t_final - t);

        if (!direct_solve)
        {
            cout << "direct == " << direct_solve << endl;
            cout << "time_dep_diff: Here?" << endl;

            ode_solver->Step(n_block, t, dt_real);
            cout << "time_dep_diff: Here?" << endl;
            ti++;

            n_e.MakeRef(&fes, n_block.GetBlock(0), 0);
            n_p.MakeRef(&fes, n_block.GetBlock(1), 0);

            // cout << "n_e check" << n_e << endl;

            done = (t >= t_final - 1e-8 * dt);

            if (done || ti % vis_steps == 0)
            {
                cout << "time step: " << ti << ", time: " << t << endl;

                dc->SetCycle(ti);
                dc->SetTime(t);
                dc->Save();
            }
        }
        else
        {
            // ode_solver->Step(n_block, t, dt_real); replace with direct

            cout << "direct == " << direct_solve << endl;

            MSp.Mult(n_block, n_temp);
            AD->Mult(n_temp, n_next);

            BlockVector n_next_block(n_next, offsets);

            n_block = n_next_block;

            ti++;

            n_e.MakeRef(&fes, n_block.GetBlock(0), 0);
            n_p.MakeRef(&fes, n_block.GetBlock(1), 0);

            // cout << "n_e check" << n_e << endl;

            done = (t >= t_final - 1e-8 * dt);

            if (done || ti % vis_steps == 0)
            {
                cout << "time step: " << ti << ", time: " << t << endl;

                dc->SetCycle(ti);
                dc->SetTime(t);
                dc->Save();
            }

            t += dt_real;
        }
    }

    n_e.MakeRef(&fes, n_block.GetBlock(0), 0);
    n_p.MakeRef(&fes, n_block.GetBlock(1), 0);

    // cout << "n_e final:"<<endl;
    // cout<< n_e << endl;
    /*
        dc = NULL;
        dc = new VisItDataCollection("Time_dep_final_1V",&mesh);
        dc->SetPrecision(8);
        dc->RegisterField("solution_e",&n_e);
        dc->RegisterField("solution_p",&n_p);
        dc->Save();
        */
    delete ode_solver;
    delete dc;
    delete AD;
}
