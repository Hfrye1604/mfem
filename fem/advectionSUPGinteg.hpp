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
 * Created on June 23, 2023, 1:00 PM
 * 
 * 
 */

#ifndef MFEM_ADVECTSUPGINTEG
#define MFEM_ADVECTSUPGINTEG

#include "../config/config.hpp"
#include "nonlininteg.hpp"
#include "fe.hpp"
#include "fespace.hpp"
#include "coefficient.hpp"
#include "ceed/interface/util.hpp"

namespace mfem
{
///Advection field integrator with SUPG
class AdvectionSUPGIntegrator : public BilinearFormIntegrator
{
private:    
    DenseMatrix dshape, adjJ,bdfVec_ir; //PMatI, PMatO
    
    VectorCoefficient* bdfVec;
  
    Coefficient* TFact;
    Coefficient* DiffCoef; //diffusion coefficient
    
    double eleLength , eleVol;
    
    Vector shape, vecPhys, BdFidxT ,BdFidxT1, BdFidxT2;
/*
    void StabConvGradIntegrator(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    DenseMatrix &elmat);
   */

    void ConvectionIntegrator(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    DenseMatrix &elmat);

    void StabSUPGIntegrator(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    DenseMatrix &elmat);
   
    void CalculateTau(const double normVel, const double Diff, double& tauSUPG);
    
public:
    AdvectionSUPGIntegrator(VectorCoefficient &advectionVelocity, 
          double TauFactor, ConstantCoefficient diffusion) : bdfVec(&advectionVelocity){
           TFact = new ConstantCoefficient(TauFactor);
           DiffCoef = new ConstantCoefficient(diffusion);
          }
          
    virtual void AssembleElementMatrix(const FiniteElement &el,
                                      ElementTransformation &Tr,
                                      DenseMatrix &elmat); //.const Vector &elfun
    
    virtual ~AdvectionSUPGIntegrator();
   
};
}


#endif /* ADVECTIONSUPGINTEG_HPP */

