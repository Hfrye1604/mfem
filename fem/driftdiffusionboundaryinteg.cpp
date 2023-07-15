#include "fem.hpp"
#include "cmath"

namespace mfem
{
    void VectorNormedMassIntegrator::AssembleFaceMatrix(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Trans,
                                   DenseMatrix &elmat)
        {
            MFEM_ASSERT(Trans.Elem2No < 0,
               "support for interior faces is not implemented");

            int nd1 = el1.GetDof();
            int dim = el1.GetDim();
            double w;

            Vector normal(dim);
            Vector vecRef(dim);
            Vector vecPhys(dim);

            #ifdef MFEM_THREAD_SAFE
            Vector shape;
            #endif

            elmat.SetSize(nd1);
            shape.SetSize(nd1);

            const IntegrationRule *ir = IntRule;
            if (ir == NULL)
            {
                int order = 2 * el1.GetOrder();

                ir = &IntRules.Get(Trans.GetGeometryType(), order);
            }

            Q->Eval(Q_ir,Trans, *ir);

            elmat = 0.0;
            for(int i = 0; i < ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);

                // Set the integration point in the face and the neighboring element
                Trans.SetAllIntPoints(&ip);

                // Access the neighboring element's integration point
                const IntegrationPoint &eip = Trans.GetElement1IntPoint();
                el1.CalcShape(eip, shape);

                CalcOrtho(Trans.Jacobian(), normal); //may need to apply jacobian
                normal *= ip.weight;

                CalcAdjugate(Trans.Jacobian(), adjJ);
                Q_ir.GetColumnReference(i, vecRef);

                adjJ.Mult(vecRef, vecPhys);

                w = vecPhys*normal; 

                double k = Factor->Eval(Trans,ip);
                w *= k;
/*
                w = Trans.Weight() * ip.weight;
                if (Q)
                {
                    w *= Q -> Eval(Trans, ip);
                }
*/
                AddMult_a_VVt(w, shape, elmat);
            }

        };

        VectorNormedMassIntegrator::~VectorNormedMassIntegrator(){}
}
