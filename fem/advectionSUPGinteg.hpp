// Copyright (c) 2023, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

/* 
 * File:   advectionSUPGinteg.hpp
 * Author: Hfrye1604
 *
 * Created on September 13, 2017, 10:55 AM
 */

#ifndef MFEM_ADVECTSUPGINTEG
#define MFEM_ADVECTSUPGINTEG

#include "../config/config.hpp"
#include "fe.hpp"
#include "coefficient.hpp"
#include "nonlininteg.hpp"

namespace mfem
{
///Navier-Stokes field integrator 
class VectorAdvectionIntegrator : public BilinearFormIntegrator
{
private:
    DenseMatrix PMatI, PMatO, dshape, gshape, Jinv;
    
    VectorCoefficient* bdfVec;
  
    Coefficient* nuCoef; 
    
    double eleVol, eleLength;
    
    Vector shape;
    
    void ConvectionIntegrator(const FiniteElement &el,
                                    ElementTransformation &Trans,
                                    DenseMatrix& elmat);
    
    void StabLaplacianIntegrator(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      DenseMatrix &elmat);
    void SatbDivIntegrator(const FiniteElement &el, 
                           ElementTransformation &Tr,
                           DenseMatrix &elmat);
    
    void StabConvectionIntegrator(const FiniteElement &el,
                                  ElementTransformation &Tr, 
                                  DenseMatrix &elmat);
    
    void StabConvGradIntegrator(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    DenseMatrix &elmat);
    
    void StabVecGradIntegrator(const FiniteElement &el,
                                       ElementTransformation &Tr,
                                       Vector &elvect);
    
    void StabVecAdvIntegrator(const FiniteElement &el,
                                      ElementTransformation &Tr, 
                                      Vector &elvect);
    
    void StabVecLapIntegrator(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      Vector &elvect);
    
    void CalcPhysLaplacian(DenseMatrix &Hessian,  
                           double nnodes,
                           double spaceDim,
                           Vector& Laplacian);
   
    void StabLapLapIntegrator(const FiniteElement &el,
                              ElementTransformation &Tr,
                              DenseMatrix &elmat);
    
    void StabLapGradIntegrator(const FiniteElement &el, 
                               ElementTransformation &Tr, 
                               DenseMatrix &elmat); 
    
    void StabLapConvIntegrator(const FiniteElement &el,
                               ElementTransformation &Tr, 
                               DenseMatrix &elmat);
    
    void CalculateTaus(const double nu, const double normVel, double& tauMom, 
                       double& tauMass);
    
public:
    VectorAdvectionIntegrator(Vector exBodyForce, 
          double ViscCoef ) {
        bdfVec = new VectorConstantCoefficient(exBodyForce);
        nuCoef = new ConstantCoefficient(ViscCoef);
          };
          
    virtual void AssembleElementVector(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      const Vector &elfun, Vector &elvect);

    virtual void AssembleElementGrad(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    const Vector &elfun, DenseMatrix &elmat);

   virtual ~VectorAdvectionIntegrator();
   
};
}


#endif /* ADVECTIONSUPGINTEG_HPP */

