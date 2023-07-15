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
        const FiniteElement &el, ElementTransformation &Trans,// const Vector &elfun,
        DenseMatrix &elmat)
        {
           int nd = el.GetDof();
            int dim = el.GetDim();

            elmat.SetSize(nd);
            dshape.SetSize(nd,dim);
            adjJ.SetSize(dim);
            shape.SetSize(nd);
            vecPhys.SetSize(dim);
            BdFidxT.SetSize(nd);
            BdFidxT1.SetSize(nd);
            BdFidxT2.SetSize(nd);

            elmat = 0.0;

            //MFEM master branch has integrator that computes the advection term
            //ConservativeConvectionIntegrator -alpha (u, q . grad v), negative transpose of ConvectionIntegrator
            //convec_term.AddDomainIntegrator(new ConservativeConvectionIntegrator(*bdfVec,-1.0));
            ConvectionIntegrator(el,Trans,elmat);

            //Compute element parameters used in Tau's calculations 
            eleVol = Geometry::Volume[el.GetGeomType()] * Trans.Weight();
            eleLength = ((dim == 3) ? (0.60046878 * pow(eleVol,0.333333333333333333333))
                : (1.128379167 * sqrt(eleVol)));    

            StabSUPGIntegrator(el,Trans,elmat);

        }

/*
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
*/

        void AdvectionSUPGIntegrator::ConvectionIntegrator(
            const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat)
        {
            Vector vecRef;

            const IntegrationRule *ir = IntRule;
            if (ir == NULL)
            {
                int order = Trans.OrderGrad(&el) + Trans.Order() + el.GetOrder();
                ir = &IntRules.Get(el.GetGeomType(), order);
            }

            bdfVec->Eval(bdfVec_ir,Trans, *ir);

            //elmat = 0.0;
            for (int i = 0; i < ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);
                el.CalcDShape(ip, dshape);
                el.CalcShape(ip, shape);

                Trans.SetIntPoint(&ip);
                CalcAdjugate(Trans.Jacobian(), adjJ);
                bdfVec_ir.GetColumnReference(i, vecRef);
                //double nu = nuCoef->Eval(Trans,ip);
                vecRef *= ip.weight;

                adjJ.Mult(vecRef, vecPhys);
                dshape.Mult(vecPhys, BdFidxT);

                AddMultVWt(shape, BdFidxT, elmat);
            }
        }

        void AdvectionSUPGIntegrator::StabSUPGIntegrator(const FiniteElement &el,
                                    ElementTransformation &Trans,
                                    DenseMatrix &elmat)
        {/*
            int nd = el.GetDof();
            int dim = el.GetDim();

            elmat.SetSize(nd);
            dshape.SetSize(nd,dim);
            adjJ.SetSize(dim);
            shape.SetSize(nd);
            vecPhys.SetSize(dim);
            BdFidxT1.SetSize(nd);
            BdFidxT1.SetSize(nd);
*/
            Vector vecRef;

            const IntegrationRule *ir = IntRule;

            if (ir == NULL)
            {
                int order = Trans.OrderGrad(&el) + Trans.Order() + el.GetOrder();
                ir = &IntRules.Get(el.GetGeomType(), order);
            }

            bdfVec->Eval(bdfVec_ir,Trans, *ir);

            //elmat = 0.0; //might need to change if including convection term in this class

            for(int i = 0; i < ir->GetNPoints();i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);
                el.CalcDShape(ip, dshape);
                el.CalcShape(ip, shape);

                Trans.SetIntPoint(&ip);
                CalcAdjugate(Trans.Jacobian(), adjJ);
                bdfVec_ir.GetColumnReference(i,vecRef);
                double tfact = TFact->Eval(Trans,ip);//changed from viscosity to scale factor; different purpose from original NS code
                double Diff = DiffCoef->Eval(Trans,ip);
                double vNorm = vecRef.Norml2();
                double tauSUPG;

                CalculateTau(vNorm, Diff, tauSUPG);

                dshape.Mult(vecRef,BdFidxT1);
                
                //vecRef *= tfact * ip.weight;
                adjJ.Mult(vecRef,vecPhys); //retrieves physical vector from reference space; vec1 is also normalized, scaled, and stabilized
                vecPhys *= tfact * tauSUPG * ip.weight ;
                //computing both Div terms
                dshape.Mult(vecPhys, BdFidxT1);
                dshape.Mult(vecRef, BdFidxT2);

                AddMultVWt(BdFidxT2,BdFidxT1, elmat);
            }
        }

        void AdvectionSUPGIntegrator::CalculateTau(
        const double normVel, const double Diff ,double& tauSUPG)
        {
            
            tauSUPG = 0.0;
            
            if(normVel != 0){
                double alpha = (normVel*eleLength)/(2.0*Diff);

                double epsilon = (1.0/tanh(alpha)) - (1.0/alpha);//Shibata paper
                
                tauSUPG = eleLength/ (2.0 * normVel) * epsilon ;//+ 4.0 * nu / (eleLength * eleLength);     
            }  
        }

        AdvectionSUPGIntegrator::~AdvectionSUPGIntegrator(){
        //fill in later    
        }

}
