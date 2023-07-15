/* 
 * File:   driftdiffusionboundaryinteg.hpp
 * Author: Hfrye1604
 *
 * Created on July 14, 2023, 1:00 PM
 * 
 * 
 */

#include "../config/config.hpp"
#include "bilininteg.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "coefficient.hpp"
#include "ceed/interface/util.hpp"

namespace mfem
{

class VectorNormedMassIntegrator : public MassIntegrator
{
protected:

#ifndef MFEM_THREAD_SAFE
   Vector shape, te_shape;
#endif
    DenseMatrix Q_ir, adjJ;
   VectorCoefficient *Q;
   Coefficient *Factor;
 /*  // PA extension
   const FiniteElementSpace *fespace;
   Vector pa_data;
   const DofToQuad *maps;         ///< Not owned
   const GeometricFactors *geom;  ///< Not owned
   int dim, ne, nq, dofs1D, quad1D;
*/
public:
    VectorNormedMassIntegrator(VectorCoefficient &q, Coefficient &factor)
    : Q(&q), Factor(&factor) {}

    using BilinearFormIntegrator::AssembleFaceMatrix;

    virtual void AssembleFaceMatrix(const FiniteElement &el1,
                                   const FiniteElement &el2,
                                   FaceElementTransformations &Trans,
                                   DenseMatrix &elmat);

    virtual ~VectorNormedMassIntegrator();
};
}