#include "../../mfem.hpp"
#include "cmath"
/*
namespace mfem
{
    class J_Solver : public Solver
    {
    private:
        SparseMatrix &M, &J, A;
        GMRESSolver linear_solver;
        Array<int> &offsets;
        Solver *invJ11, *invJ22;
        BlockDiagonalPreconditioner *J_prec;
        //BlockILU prec;
        double dt;

    public:
        J_Solver(SparseMatrix& M_, SparseMatrix& J_, Array<int> &offsets_ ) //, const FiniteElementSpace &fes
        : M(M_), J(J_), dt(-1.0), offsets(offsets_) //prec(fes.GetFE(0)->GetDof(),BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
             {
                int block_size = J.Size()/2; 
                DenseMatrix J11(block_size);
                
                DenseMatrix J22(block_size);

                Array<int> row11(block_size);
                Array<int> row22(block_size);
                Array<int> col11(block_size);
                Array<int> col22(block_size);

                for(int i=offsets[0]; i < offsets[1]; i++)
                    {
                        row11[i] = i;
                        col11[i] = i;
                        row22[i] = i + offsets[1];
                        col22[i] = i + offsets[1];
                    }   
           
                
                J.GetSubMatrix(row11,col11,J11);
                J.GetSubMatrix(row22,col22,J22);
std::cout <<"Here?" << std::endl;
                invJ11 = new DSmoother((SparseMatrix&) J11);
                invJ22 = new DSmoother((SparseMatrix&) J22);
std::cout <<"Here?" << std::endl;
                J_prec = new BlockDiagonalPreconditioner(offsets);
std::cout <<"Here?" << std::endl;
                J_prec->SetDiagonalBlock(0,invJ11);
                J_prec->SetDiagonalBlock(1,invJ22);
std::cout <<"Here?" << std::endl;
                linear_solver.iterative_mode = false;
                linear_solver.SetRelTol(1e-9);
                linear_solver.SetAbsTol(0.0);
                linear_solver.SetMaxIter(500);
                linear_solver.SetPrintLevel(0);
                linear_solver.SetPreconditioner(*J_prec);
std::cout <<"Here?" << std::endl;
             }

        void SetTimeStep(double dt_)
        {
            if (dt_ != dt)
            {
                dt = dt_;
                // Form operator A = M - dt*K
                A = J;
                A *= -dt;
                A += M;

                // this will also call SetOperator on the preconditioner
                linear_solver.SetOperator(A);
            }
        };

        void SetOperator(const Operator &op)
        {
            linear_solver.SetOperator(op);
        };

        virtual void Mult(const Vector &x, Vector &y) const
        {
            linear_solver.Mult(x,y);
        };

        virtual ~J_Solver(){
            //delete linear_solver;
        };
    };

    class FE_Evolution : public TimeDependentOperator
    {
    private:
        SparseMatrix &M, &J;
        //BlockMatrix &J;
        Array<int> &offsets;
        const Vector &b;
        Solver *M_prec;
        CGSolver M_solver;
        J_Solver *j_solver;

        mutable Vector z;
    public:
        FE_Evolution(SparseMatrix &M_, SparseMatrix &J_, const Vector &b_, Array<int>& offsets_)
            : TimeDependentOperator(M_.Height()), M(M_), J(J_), b(b_), z(M_.Height()), offsets(offsets_)
            {
                //Array<int> ess_tdof_list;
std::cout <<"FE evo: Here?" << std::endl;
                M_prec = new DSmoother(M);
                M_solver.SetOperator(M);
                j_solver = new J_Solver(M, J, offsets); //, *M.FESpace()
std::cout <<"Fe evo: Here?" << std::endl;
                M_solver.SetPreconditioner(*M_prec);
                M_solver.iterative_mode = false;
                M_solver.SetRelTol(1e-9);
                M_solver.SetAbsTol(0.0);
                M_solver.SetMaxIter(100);
                M_solver.SetPrintLevel(0);
std::cout <<"Fe evo: Here?" << std::endl;
            };

        virtual void Mult(const Vector &x, Vector &y) const
        {
            // y = M^{-1} (J x + b)
std::cout <<"Fe evo; Mult: Here?" << std::endl;
            J.Mult(x, z);
            z += b;
            M_solver.Mult(z,y);
std::cout <<"Fe evo; Mult: Here?" << std::endl;
        };

        virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k)
        {
std::cout <<"Fe evo; ImplicitSolve: Here?" << std::endl;
            J.Mult(x,z);
            z += b;
            j_solver->SetTimeStep(dt);
            j_solver->Mult(z,k);
std::cout <<"Fe evo; ImplicitSolve: Here?" << std::endl;
        };

        virtual ~FE_Evolution()
        {
            delete M_prec;
            delete j_solver;
        };

    };
}
*/

namespace mfem
{
    class J_Solver : public Solver
    {
    private:
        SparseMatrix &M, &J, A;
        GMRESSolver linear_solver;
        Array<int> &offsets;
        BlockILU *prec;
        double dt;
    public:
        J_Solver(SparseMatrix& M_, SparseMatrix& J_, Array<int> &offsets_)
        : M(M_), J(J_), dt(-1.0), offsets(offsets_)
             {
                prec = new BlockILU(offsets[1],BlockILU::Reordering::MINIMUM_DISCARDED_FILL);
                linear_solver.iterative_mode = false;
                linear_solver.SetRelTol(1e-9);
                linear_solver.SetAbsTol(0.0);
                linear_solver.SetMaxIter(100);
                linear_solver.SetPrintLevel(0);
                linear_solver.SetPreconditioner(*prec);
             }

        void SetTimeStep(double dt_)
        {
            if (dt_ != dt)
            {
                dt = dt_;
                // Form operator A = M - dt*K
                A = J;
                A *= -dt;
                A += M;

                // this will also call SetOperator on the preconditioner
                linear_solver.SetOperator(A);
            }
        };

        void SetOperator(const Operator &op)
        {
            linear_solver.SetOperator(op);
        };

        virtual void Mult(const Vector &x, Vector &y) const
        {
            linear_solver.Mult(x,y);
        };

        virtual ~J_Solver(){
            delete prec;
        };
    };

    class FE_Evolution : public TimeDependentOperator
    {
    private:
        SparseMatrix &M, &J;
        const Vector &b;
        Array<int> &offsets;
        Solver *M_prec;
        CGSolver M_solver;
        J_Solver *j_solver;

        mutable Vector z;
    public:
        FE_Evolution(SparseMatrix &M_, SparseMatrix &J_, const Vector &b_, Array<int> &offsets_)
            : TimeDependentOperator(M_.Height()), M(M_), J(J_), b(b_), z(M_.Height()), offsets(offsets_)
            {
                Array<int> ess_tdof_list;
                M_prec = new DSmoother(M);
                M_solver.SetOperator(M);
                j_solver = new J_Solver(M, J, offsets); //, *M.FESpace()

                M_solver.SetPreconditioner(*M_prec);
                M_solver.iterative_mode = false;
                M_solver.SetRelTol(1e-9);
                M_solver.SetAbsTol(0.0);
                M_solver.SetMaxIter(100);
                M_solver.SetPrintLevel(0);
            };

        virtual void Mult(const Vector &x, Vector &y) const
        {
            // y = M^{-1} (J x + b)
            J.Mult(x, z);
            z += b;
            M_solver.Mult(z,y);
        };

        virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k)
        {
            J.Mult(x,z);
            z += b;
            j_solver->SetTimeStep(dt);
            j_solver->Mult(z,k);
        };

        virtual ~FE_Evolution()
        {
            delete M_prec;
            delete j_solver;
        };

    };
}

/*
#include "../../mfem.hpp"
#include "cmath"

namespace mfem
{
    class J_Solver : public Solver
    {
    private:
        SparseMatrix &M, &J, A;
        GMRESSolver linear_solver;
        //BlockILU prec;
        double dt;
    public:
        J_Solver(SparseMatrix& M_, SparseMatrix& J_) //, const FiniteElementSpace &fes
        : M(M_), J(J_), dt(-1.0) //prec(fes.GetFE(0)->GetDof(),BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
             {
                linear_solver.iterative_mode = false;
                linear_solver.SetRelTol(1e-9);
                linear_solver.SetAbsTol(0.0);
                linear_solver.SetMaxIter(100);
                linear_solver.SetPrintLevel(0);
                //linear_solver.SetPreconditioner(prec);
             }

        void SetTimeStep(double dt_)
        {
            if (dt_ != dt)
            {
                dt = dt_;
                // Form operator A = M - dt*K
                A = J;
                A *= -dt;
                A += M;

                // this will also call SetOperator on the preconditioner
                linear_solver.SetOperator(A);
            }
        };

        void SetOperator(const Operator &op)
        {
            linear_solver.SetOperator(op);
        };

         void Mult(const Vector &x, Vector &y) const
        {
            linear_solver.Mult(x,y);
        };

         ~J_Solver(){
            //delete linear_solver;
        };
    };

    class FE_Evolution : public TimeDependentOperator
    {
    private:
        SparseMatrix &M, &J;
        const Vector &b;
        Solver *M_prec;
        CGSolver M_solver;
        J_Solver *j_solver;

        mutable Vector z;
    public:
        FE_Evolution(SparseMatrix &M_, SparseMatrix &J_, const Vector &b_)
            : TimeDependentOperator(M_.Height()), M(M_), J(J_), b(b_), z(M_.Height())
            {
                //z = new Vector(M.Height());
                Array<int> ess_tdof_list;
                M_prec = new DSmoother(M);
                M_solver.SetOperator(M);
                j_solver = new J_Solver(M, J); //, *M.FESpace()

                M_solver.SetPreconditioner(*M_prec);
                M_solver.iterative_mode = false;
                M_solver.SetRelTol(1e-9);
                M_solver.SetAbsTol(0.0);
                M_solver.SetMaxIter(100);
                M_solver.SetPrintLevel(0);
            };

        void Mult(const Vector &x, Vector &y) const
        {
            // y = M^{-1} (J x + b)
    
            J.Mult(x, z);
            z += b;
            M_solver.Mult(z,y);
        };

        void ImplicitSolve(const double dt, const Vector &x, Vector &k)
        {
            std::cout <<"z size: " <<z.Size() << std::endl;
            J.Mult(x,z);
            z += b;
            j_solver->SetTimeStep(dt);
            j_solver->Mult(z,k);
        };

        ~FE_Evolution()
        {
            delete M_prec;
            delete j_solver;
        };

    };
}
*/