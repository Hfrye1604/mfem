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

#include "fem.hpp"
#include "cmath"

namespace mfem
{
    void AdvectionSUPGIntegrator::AssembleElementMatrix(
        const FiniteElement &el, ElementTransformation &Tr,// const Vector &elfun,
        DenseMatrix &elmat)
        {
            int dim = el.GetDim() + 1; //1 for pressure 
            int nd  = el.GetDof();
            int sp_dim = Tr.GetSpaceDim();          
            
            dshape.SetSize(nd, sp_dim); 
            gshape.SetSize(nd, sp_dim);   
            Jinv  .SetSize(sp_dim);

//Will have to figure out alternative
            //PMatI.UseExternalData(elfun.GetData(), nd, dim);
            elmat.SetSize(nd*dim);
            //PMatO.UseExternalData(elvect.GetData(), nd, dim);
            
            elmat = 0.0;

            //MFEM master branch has integrator that computes the advection term
            //ConservativeConvectionIntegrator -alpha (u, q . grad v), negative transpose of ConvectionIntegrator
            //convec_term.AddDomainIntegrator(new ConservativeConvectionIntegrator(*bdfVec,-1.0));

            //Compute element parameters used in Tau's calculations 
            eleVol = Geometry::Volume[el.GetGeomType()] * Tr.Weight();
            eleLength = ((sp_dim == 3) ? (0.60046878 * pow(eleVol,0.333333333333333333333))
                : (1.128379167 * sqrt(eleVol)));    

            //StabLaplacianIntegrator(el,Tr,elmat);
            StabConvGradIntegrator(el,Tr,elmat);

        }


        void AdvectionSUPGIntegrator::StabConvGradIntegrator(const FiniteElement &el,
                                    ElementTransformation &Tr,
                                    DenseMatrix &elmat)
        {
            int dim = el.GetDim() + 1;
            int nd  = el.GetDof();
            int sp_dim = Tr.GetSpaceDim();
            DenseMatrix pelemat;
            Vector AdvUP(dim), advGrad(nd), auxU(sp_dim), vec1;
            double norm;
            
            pelemat.SetSize(nd);
            
            int intorder = el.GetOrder() + 2 * Tr.OrderGrad(&el);
            const IntegrationRule &ir = IntRules.Get(el.GetGeomType(), intorder);
            
            for (int i = 0; i < ir.GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir.IntPoint(i);
                
                el.CalcDShape(ip, dshape);
                el.CalcShape (ip, shape);
                
                Tr.SetIntPoint (&ip);
                norm = ip.weight * Tr.Weight();
                CalcInverse(Tr.Jacobian(), Jinv);
                
                Mult (dshape, Jinv, gshape);
                
                //compute tau
                double tauSUPG;
                double nu = nuCoef->Eval (Tr, ip);
                double Diff = DiffCoef->Eval (Tr,ip);

                PMatI.MultTranspose(shape, AdvUP);
                for (int kk = 0; kk < sp_dim; kk++)
                auxU[kk] = AdvUP[kk]; 
            
                double Unorm = auxU.Norml2();
                CalculateTaus(nu, Unorm, Diff ,tauSUPG);      
                norm *= tauSUPG;   
                
                gshape.Mult(auxU, advGrad);

                for (int kk = 0; kk < sp_dim; kk++)
                {
                    gshape.GetColumnReference(kk, vec1);
                    MultVWt(advGrad, vec1, pelemat);
                    pelemat *= norm;
                    elmat.AddMatrix(pelemat, nd * kk, sp_dim * nd);
                    pelemat.Transpose();
                    elmat.AddMatrix(pelemat, sp_dim * nd, nd * kk);
                }         
            }
        }   

        void AdvectionSUPGIntegrator::CalculateTaus(const double nu, 
        const double normVel, const double Diff ,double& tauSUPG)
        {
            
            tauSUPG = 0.0;

            double alpha = (normVel*eleLength)/(2.0*Diff);

            double epsilon = (1.0/tanh(alpha)) - (1.0/alpha);//Shibata paper
            
            double tauSUPG = eleLength/ (2.0 * normVel) * epsilon ;//+ 4.0 * nu / (eleLength * eleLength);     
               
        }

        ~AdvectionSUPGIntegrator(){
        //fill in later    
        }

}
