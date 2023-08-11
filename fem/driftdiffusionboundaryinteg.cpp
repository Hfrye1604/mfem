#include "fem.hpp"
#include "cmath"

namespace mfem
{
    void VectorNormedMassIntegrator::AssembleElementMatrix(const FiniteElement &el,
                                   ElementTransformation &Trans,
                                   DenseMatrix &elmat)
        {
            //MFEM_ASSERT(Trans.Elem2No < 0,
             //  "support for interior faces is not implemented");

            int nd = el.GetDof();
            int dim = Trans.GetSpaceDim();
            double w;

            //std::cout << "fem space dimension: "<<dim << std::endl;

            Vector normal(dim);
            Vector vecRef(dim);
            Vector vecPhys(dim);

            #ifdef MFEM_THREAD_SAFE
            Vector shape;
            #endif

            elmat.SetSize(nd);
            shape.SetSize(nd);

            const IntegrationRule *ir = IntRule;
            if (ir == NULL)
            {
                int order = Trans.OrderGrad(&el) + Trans.Order() + el.GetOrder();
                ir = &IntRules.Get(el.GetGeomType(), order);
            }
            Q->Eval(Q_ir,Trans, *ir);

            elmat = 0.0;
            for(int i = 0; i < ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);

                Trans.SetIntPoint(&ip);
                
                el.CalcPhysShape(Trans, shape);

                if(dim > 1){
                    CalcOrtho(Trans.Jacobian(), normal); //may need to apply jacobian
                    normal *= ip.weight;
                    Q_ir.GetColumnReference(i, vecRef);
                    w = vecRef*normal; 
                    double k = Factor->Eval(Trans,ip);
                    w *= k;
                    AddMult_a_VVt(w, shape, elmat);
                } else{
                    normal *= ip.weight;
                    Q_ir.GetColumnReference(i, vecRef);
                    w = vecRef*normal; 
                    double k = Factor->Eval(Trans,ip);
                    w *= k;
                    AddMult_a_VVt(w, shape, elmat);
                    normal *= -1;
                }

            }

        };

        VectorNormedMassIntegrator::~VectorNormedMassIntegrator(){};
}
