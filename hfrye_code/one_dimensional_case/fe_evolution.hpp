#include "mfem.hpp"
#include "cmath"

namespace mfem
{
    class J_Solver : public Solver
    {
    private:
        SparseMatrix &M, &J, A;
        GMRESSolver linear_solver;
        BlockILU prec;
        double dt;
    public:
        J_Solver(SparseMatrix& M_, SparseMatrix& J_, const FiniteElementSpace &fes)
        : M(M_), J(J_), prec(fes.GetFE(0)->GetDof(), //might have to double check to see if this is what we want
             BlockILU::Reordering::MINIMUM_DISCARDED_FILL), dt(-1.0) 
             {
                linear_solver.iterative_mode = false;
                linear_solver.SetRelTol(1e-9);
                linear_solver.SetAbsTol(0.0);
                linear_solver.SetMaxIter(100);
                linear_solver.SetPrintLevel(0);
                linear_solver.SetPreconditioner(prec);
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
        }

        void SetOperator(const Operator &op)
        {
            linear_solver.SetOperator(op);
        }

        virtual void Mult(const Vector &x, Vector &y) const
        {
            linear_solver.Mult(x,y);
        }
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
                Array<int> ess_tdof_list;
                M_prec = new DSmoother(M);
                M_solver.SetOperator(M);
                j_solver = new J_Solver(M, J, *M.FESpace());

                M_solver.SetPreconditioner(*M_prec);
                M_solver.iterative_mode = false;
                M_solver.SetRelTol(1e-9);
                M_solver.SetAbsTol(0.0);
                M_solver.SetMaxIter(100);
                M_solver.SetPrintLevel(0);

            };

        virtual void Mult(const Vector &x, Vector &y) const
        {
            // y = M^{-1} (K x + b)
            J.Mult(x, z);
            z += b;//#include <lapack>
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