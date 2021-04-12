
#include "integrator.hpp"

#ifdef MFEM_USE_MPI

void
EHSMaxwellSolver::Init(mfem::Operator &P, mfem::TimeDependentOperator & F)
{
   P_ = &P; F_ = &F;

   dE_.SetSize(F_->Height());
   dB_.SetSize(P_->Height());
}

EHSVMaxwellSolver::EHSVMaxwellSolver(int order)
   : order_(order)
{
   a_.SetSize(order);
   b_.SetSize(order);

   switch (order_)
   {
      case 1:
         a_[0] = 1.0;
         b_[0] = 1.0;
         break;
      case 2:
         a_[0] = 0.5;
         a_[1] = 0.5;
         b_[0] = 0.0;
         b_[1] = 1.0;
         break;
      case 4:
         a_[0] = 0.5/(2.0-std::pow(2.0,1.0/3.0));
         a_[1] = 0.5 - a_[0];
         a_[2] = a_[1];
         a_[3] = a_[0];
         b_[0] = 0.0;
         b_[1] = 1.0/(2.0-std::pow(2.0,1.0/3.0));
         b_[2] = 1.0 - 2.0 * b_[1];
         b_[3] = b_[1];
         break;
      default:
         MFEM_ASSERT(false, "Unsupported order in EHSVSolver");
   };
}

void
EHSVMaxwellSolver::Step(mfem::Vector &B, mfem::Vector &E, double &t, double &dt)
{
   for (int i=0; i<order_; i++)
   {
      if ( b_[i] != 0.0 )
      {
         F_->SetTime(t);
         if ( F_->isExplicit() )
         {
            F_->Mult(B, dE_);
         }
         else
         {
            F_->ImplicitSolve(b_[i] * dt, B, dE_);
         }
         E.Add(b_[i] * dt, dE_);
      }

      P_->Mult(E, dB_);
      B.Add(a_[i] * dt, dB_);

      t += a_[i] * dt;
   }
}

EHSVSolver::EHSVSolver(int order, VMSolver &vmsolver, mfem::Operator & P, mfem::TimeDependentOperator & F)
   : order_(order),
     vmsolver_(&vmsolver),
     F_(&vmsolver),
     P_(&vmsolver.GetNegCurl())
{
   dE_.SetSize(F_->Height());
   dB_.SetSize(P_->Height());
   switch (order_)
   {
      case 1:
         blocks_ = 1;
         a_.SetSize(blocks_*4);
         a_[0] = 1.0;
         a_[1] = 1.0;
         a_[2] = 1.0;
         a_[3] = 0.0;
         a_[4] = 0.0;
         break;
      case 2:
         blocks_ = 1;
         a_.SetSize(blocks_*4);
         a_[0] = 0.0;
         a_[1] = 0.0;
         a_[2] = 0.5;
         a_[3] = 1.0;
         break;
      case 4:
         blocks_ = 3;
         a_.SetSize(blocks_*4);
         alpha_ = 1.0/(2.0-std::pow(2.0,1.0/3.0));
         beta_ = 1.0 - 2.0*alpha_;
         a_[0] = 0.0*alpha_;
         a_[1] = 0.0*alpha_;
         a_[2] = 0.5*alpha_;
         a_[3] = 1.0*alpha_;
         a_[4] = 0.0*beta_;
         a_[5] = 0.0*beta_;
         a_[6] = 0.5*beta_;
         a_[7] = 1.0*beta_;
         a_[8] = 0.0*alpha_;
         a_[9] = 0.0*alpha_;
         a_[10] = 0.5*alpha_;
         a_[11] = 1.0*alpha_;
         break;
      default:
         MFEM_ASSERT(false, "Unsupported order in EHSVSolver");
   };
}

void
EHSVSolver::Step1(mfem::Vector &B, mfem::Vector &E, double &t, double &dt)
{
   if ( a_[0] != 0.0 )
   {
      F_->SetTime(t);
      if ( F_->isExplicit() )
      {
         F_->Mult(B, dE_);
      }
      else
      {
         F_->ImplicitSolve(a_[0] * dt, B, dE_);
      }
      E.Add(a_[0] * dt, dE_);
   }
   if ( a_[1] != 0.0 )
   {
      P_->Mult(E, dB_);
      B.Add(a_[1] * dt, dB_);
      t += a_[1] * dt;
   }
   vmsolver_->fields_->SyncGridFuncs();
   if ( a_[2] != 0.0 )
   {
      Push_x1(a_[2]*dt);
   }
   /*if ( a_[3] != 0.0 )
   {
      //push y
      //push z
   }*/
}

void
EHSVSolver::Step2(mfem::Vector &B, mfem::Vector &E, double &t, double &dt)
{

   /*if ( a_[0] != 0.0 )
   {
      //push x
   }
   if ( a_[1] != 0.0 )
   {
      //push y
      //push z
   }*/
   if ( a_[2] != 0.0 )
   {
      P_->Mult(E, dB_);
      B.Add(a_[2] * dt, dB_);
      t += a_[2] * dt;
   }
   if ( a_[3] != 0.0 )
   {
      F_->SetTime(t);
      if ( F_->isExplicit() )
      {
         F_->Mult(B, dE_);
      }
      else
      {
         F_->ImplicitSolve(a_[3] * dt, B, dE_);
      }
      E.Add(a_[3] * dt, dE_);
   }
   if ( a_[2] != 0.0 )
   {
      P_->Mult(E, dB_);
      B.Add(a_[2] * dt, dB_);
      t += a_[2] * dt;
   }
   vmsolver_->fields_->SyncGridFuncs();
   /*if ( a_[1] != 0.0 )
   {
      //push z
      //push y
   }*/
}


void
EHSVSolver::Step4(mfem::Vector &B, mfem::Vector &E, double &t, double &dt)
{
   for (int i=0; i<blocks_; i++)
   {
      /*if ( a_[4*i] != 0.0 )
      {
         //push x
      }
      if ( a_[4*i+1] != 0.0 )
      {
         //push y
         //push z
      }*/
      if ( a_[4*i+2] != 0.0 )
      {
         P_->Mult(E, dB_);
         B.Add(a_[4*i+2] * dt, dB_);
         t += a_[4*i+2] * dt;
      }
      if ( a_[4*i+3] != 0.0 )
      {
         F_->SetTime(t);
         if ( F_->isExplicit() )
         {
            F_->Mult(B, dE_);
         }
         else
         {
            F_->ImplicitSolve(a_[4*i+3] * dt, B, dE_);
         }
         E.Add(a_[4*i+3] * dt, dE_);
      }
      if ( a_[4*i+2] != 0.0 )
      {
         P_->Mult(E, dB_);
         B.Add(a_[4*i+2] * dt, dB_);
         t += a_[4*i+2] * dt;
      }
      vmsolver_->fields_->SyncGridFuncs();
      /*if ( a_[4*i+1] != 0.0 )
      {
         //push z
         //push y
      }*/
   }
}

void
EHSVSolver::Run(mfem::Vector &B, mfem::Vector &E, double &t, double &dt, double tf)
{
   switch (order_)
   {
      case 1:
         while (t < tf) { Step1(B, E, t, dt); }
         break;
      case 2:
         Step2(B, E, t, dt);
         a_[0] *= 2.0;
         while (t < tf) { Step2(B, E, t, dt); }
         a_[0] *= 0.5;
         //push x
         break;
      case 4:
         Step4(B, E, t, dt);
         a_[0] *= 2.0;
         while (t < tf) { Step4(B, E, t, dt); }
         a_[0] *= 0.5;
         //push x
         break;
      default:
         MFEM_ASSERT(false, "Unsupported order in EHSVSolver");
   };
   //initial full step with first alpha/2
   //steps with first alpha
   //final thetax with alpha/2
}

void
EHSVSolver::Push_x1(double dt)
{
   double dx1, dE1;

   for (int i=0; i< vmsolver_->fields_->pmesh_->GetNE(); i++)
   {
      for (int j=0; j<vmsolver_->param_->num_species; j++)
      {

         VMParticles* partsInElem = (*vmsolver_->particles_->particlesArray_)(i,j);

         int npts = partsInElem->size;

         mfem::ElementTransformation * etrans = vmsolver_->fields_->H1FESpace_->GetElementTransformation(i);

         for (int k=0; k<npts; k++)
         {
            mfem::IntegrationPoint ip;
            mfem::Vector pvec;
            ip.Set3((*partsInElem->x1_)(k),(*partsInElem->x2_)(k),(*partsInElem->x3_)(k));
            etrans->SetIntPoint(&ip);
            mfem::DenseMatrix invJacobian = etrans->InverseJacobian();
            dx1 = vmsolver_->param_->N_v_s[j]*invJacobian(0,0)*(*partsInElem->v1_)(k)*dt;
            if ((*partsInElem->x1_)(k)+dx1<0 || (*partsInElem->x1_)(k)+dx1>1)
            {
               std::cout << "dx1: " << (*partsInElem->x1_)(k)+dx1 << std::endl;
            }
         }
      }
   }
}

void
EHSVSolver::Push_x2(double dt)
{

   for (int i=0; i< vmsolver_->fields_->pmesh_->GetNE(); i++)
   {
      for (int j=0; j<vmsolver_->param_->num_species; j++)
      {

         VMParticles* partsInElem = (*vmsolver_->particles_->particlesArray_)(i,j);

         int npts = partsInElem->size;

         mfem::ElementTransformation * etrans = vmsolver_->fields_->H1FESpace_->GetElementTransformation(i);

         for (int k=0; k<npts; k++)
         {
            mfem::IntegrationPoint ip;
            mfem::Vector pvec;
            ip.Set3((*partsInElem->x1_)(k),(*partsInElem->x2_)(k),(*partsInElem->x3_)(k));
            etrans->SetIntPoint(&ip);
            mfem::DenseMatrix invJacobian = etrans->InverseJacobian();

         }
      }
   }
}

void
EHSVSolver::Push_x3(double dt)
{

   for (int i=0; i< vmsolver_->fields_->pmesh_->GetNE(); i++)
   {
      for (int j=0; j<vmsolver_->param_->num_species; j++)
      {

         VMParticles* partsInElem = (*vmsolver_->particles_->particlesArray_)(i,j);

         int npts = partsInElem->size;

         mfem::ElementTransformation * etrans = vmsolver_->fields_->H1FESpace_->GetElementTransformation(i);

         for (int k=0; k<npts; k++)
         {
            mfem::IntegrationPoint ip;
            mfem::Vector pvec;
            ip.Set3((*partsInElem->x1_)(k),(*partsInElem->x2_)(k),(*partsInElem->x3_)(k));
            etrans->SetIntPoint(&ip);
            mfem::DenseMatrix invJacobian = etrans->InverseJacobian();
         }
      }
   }
}



#endif // MFEM_USE_MPI
