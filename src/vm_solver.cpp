
#include "vm_solver.hpp"


#ifdef MFEM_USE_MPI


// Used for combining scalar coefficients
double prodFunc(double a, double b) { return a * b; }

VMSolver::VMSolver(ParametersInput & param, VMFields & fields, VMParticlesArray & particles)
   : dtMax_(-1.0),
     dtScale_(1.0e6),
     param_(&param),
     fields_(&fields),
     particles_(&particles),
     hDivMassMuInv_(NULL),
     hCurlLosses_(NULL),
     weakCurlMuInv_(NULL),
     Curl_(NULL),
     rhs_(NULL),
     M1Losses_(NULL),
     M2MuInv_(NULL),
     NegCurl_(NULL),
     WeakCurlMuInv_(NULL),
     HD_(NULL),
     RHS_(NULL)
{

   // Require implicit handling of loss terms
   type = fields_->lossy_ ? IMPLICIT : EXPLICIT;

   this->height = fields_->HCurlFESpace_->GlobalTrueVSize();
   this->width  = fields_->HDivFESpace_->GlobalTrueVSize();

   // Bilinear Forms
   if ( fields_->myid_ == 0 && fields_->logging_ > 0 )
   {
      std::cout << "Creating H(Div) Mass Operator" << std::endl;
   }
   hDivMassMuInv_ = new mfem::ParBilinearForm(fields_->HDivFESpace_);
   hDivMassMuInv_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*fields_->curlBCoef_));

   if ( fields_->myid_ == 0 && fields_->logging_ > 0 )
   {
      std::cout << "Creating Weak Curl Operator" << std::endl;
   }
   weakCurlMuInv_ = new mfem::ParMixedBilinearForm(fields_->HDivFESpace_,fields_->HCurlFESpace_);
   weakCurlMuInv_->AddDomainIntegrator(
      new mfem::MixedVectorWeakCurlIntegrator(*fields_->curlBCoef_));

   // Assemble Matrices
   hDivMassMuInv_->Assemble();
   weakCurlMuInv_->Assemble();

   hDivMassMuInv_->Finalize();
   weakCurlMuInv_->Finalize();

   if ( fields_->sigmaCoef_ || fields_->etaInvCoef_ )
   {
      if ( fields_->myid_ == 0 && fields_->logging_ > 0 )
      {
         std::cout << "Creating H(Curl) Loss Operator" << std::endl;
      }
      hCurlLosses_  = new mfem::ParBilinearForm(fields_->HCurlFESpace_);
      if ( fields_->sigmaCoef_ )
      {
         if ( fields_->myid_ == 0 && fields_->logging_ > 0 )
         {
            std::cout << "Adding domain integrator for conductive regions" << std::endl;
         }
         hCurlLosses_->AddDomainIntegrator(
            new mfem::VectorFEMassIntegrator(*fields_->sigmaCoef_));
      }
      if ( fields_->etaInvCoef_ )
      {
         if ( fields_->myid_ == 0 && fields_->logging_ > 0 )
         {
            std::cout << "Adding boundary integrator for absorbing boundary" << std::endl;
         }
         hCurlLosses_->AddBoundaryIntegrator(
            new mfem::VectorFEMassIntegrator(*fields_->etaInvCoef_), fields_->abc_marker_);
      }
      hCurlLosses_->Assemble();
      hCurlLosses_->Finalize();
      M1Losses_ = hCurlLosses_->ParallelAssemble();
   }

   // Create Linear Algebra Matrices
   M2MuInv_ = hDivMassMuInv_->ParallelAssemble();
   WeakCurlMuInv_ = weakCurlMuInv_->ParallelAssemble();

   if ( fields_->myid_ == 0 && fields_->logging_ > 0 )
   {
      std::cout << "Creating discrete curl operator" << std::endl;
   }
   Curl_ = new mfem::common::ParDiscreteCurlOperator(fields_->HCurlFESpace_, fields_->HDivFESpace_);
   Curl_->Assemble();
   Curl_->Finalize();
   NegCurl_ = Curl_->ParallelAssemble();
   // Beware this modifies the matrix stored within the Curl_ object.
   *NegCurl_ *= -1.0*param_->N_curlE;

   if ( fields_->myid_ == 0 && fields_->logging_ > 0 )
   {
      std::cout << "Created discrete curl operator" << std::endl;
   }

   // Build grid functions
   rhs_  = new mfem::ParGridFunction(fields_->HCurlFESpace_);

   HD_  = new mfem::HypreParVector(fields_->HDivFESpace_);
   RHS_ = new mfem::HypreParVector(fields_->HCurlFESpace_);

   dtMax_ = GetMaximumTimeStep();
}

VMSolver::~VMSolver()
{
   delete HD_;
   delete RHS_;

   delete rhs_;

   delete Curl_;

   delete M1Losses_;
   delete M2MuInv_;
   delete NegCurl_;
   delete WeakCurlMuInv_;

   delete hDivMassMuInv_;
   delete hCurlLosses_;
   delete weakCurlMuInv_;

   std::map<int, mfem::ParBilinearForm*>::iterator mit1;
   for (mit1=a1_.begin(); mit1!=a1_.end(); mit1++)
   {
      int i = mit1->first;
      delete pcg_[i];
      delete diagScale_[i];
      delete A1_[i];
      delete a1_[i];
   }

   std::map<int, mfem::Coefficient*>::iterator mit2;
   for (mit2=dtCoef_.begin(); mit2!=dtCoef_.end(); mit2++)
   {
      delete mit2->second;
   }
   for (mit2=dtSigmaCoef_.begin(); mit2!=dtSigmaCoef_.end(); mit2++)
   {
      delete mit2->second;
   }
   for (mit2=dtEtaInvCoef_.begin(); mit2!=dtEtaInvCoef_.end(); mit2++)
   {
      delete mit2->second;
   }

}


void
VMSolver::Mult(const mfem::Vector &B, mfem::Vector &dEdt) const
{
   implicitSolve(0.0, B, dEdt);
}

void
VMSolver::ImplicitSolve(double dt, const mfem::Vector &B, mfem::Vector &dEdt)
{
   implicitSolve(dt, B, dEdt);
}

void
VMSolver::setupSolver(const int idt, const double dt) const
{
   if ( pcg_.find(idt) == pcg_.end() )
   {
      if ( fields_->myid_ == 0 && fields_->logging_ > 0 )
      {
         std::cout << "Creating implicit operator for dt = " << dt << std::endl;
      }

      a1_[idt] = new mfem::ParBilinearForm(fields_->HCurlFESpace_);
      a1_[idt]->AddDomainIntegrator(
         new mfem::VectorFEMassIntegrator(fields_->epsCoef_));

      if ( idt != 0 )
      {
         dtCoef_[idt] = new mfem::ConstantCoefficient(0.5 * dt);
         if ( fields_->sigmaCoef_ )
         {
            dtSigmaCoef_[idt] = new mfem::TransformedCoefficient(dtCoef_[idt],
                                                           fields_->sigmaCoef_,
                                                           prodFunc);
            a1_[idt]->AddDomainIntegrator(
               new mfem::VectorFEMassIntegrator(dtSigmaCoef_[idt]));
         }
         if ( fields_->etaInvCoef_ )
         {
            dtEtaInvCoef_[idt] = new mfem::TransformedCoefficient(dtCoef_[idt],
                                                            fields_->etaInvCoef_,
                                                            prodFunc);
            a1_[idt]->AddBoundaryIntegrator(
               new mfem::VectorFEMassIntegrator(dtEtaInvCoef_[idt]),
               const_cast<mfem::Array<int>&>(fields_->abc_marker_));
         }
      }

      a1_[idt]->Assemble();
      a1_[idt]->Finalize();
      A1_[idt] = a1_[idt]->ParallelAssemble();

      diagScale_[idt] = new mfem::HypreDiagScale(*A1_[idt]);
      pcg_[idt] = new mfem::HyprePCG(*A1_[idt]);
      pcg_[idt]->SetTol(1.0e-12);
      pcg_[idt]->SetMaxIter(200);
      pcg_[idt]->SetPrintLevel(0);
      pcg_[idt]->SetPreconditioner(*diagScale_[idt]);
   }
}

void
VMSolver::implicitSolve(double dt, const mfem::Vector &B, mfem::Vector &dEdt) const
{
   int idt = hCurlLosses_ ? ((int)(dtScale_ * dt / dtMax_)) : 0;

   fields_->b_->Distribute(B);
   weakCurlMuInv_->Mult(*fields_->b_, *rhs_);

   if ( hCurlLosses_ )
   {
      fields_->e_->Distribute(*fields_->E_);
      hCurlLosses_->AddMult(*fields_->e_, *rhs_, -1.0);
   }

   if ( fields_->jd_ )
   {
      fields_->jCoef_->SetTime(t); // 't' is member data from mfem::TimeDependentOperator
      fields_->jd_->Assemble();
      *rhs_ -= *fields_->jd_;
   }

   if ( fields_->dEdtBCCoef_ )
   {
      fields_->dEdtBCCoef_->SetTime(t);
      fields_->dedt_->ProjectBdrCoefficientTangent(*fields_->dEdtBCCoef_,
                                          const_cast<mfem::Array<int>&>(fields_->dbc_marker_));
   }

   // Create objects and matrices for solving with the given time step
   setupSolver(idt, dt);

   // Apply essential BCs and determine true DoFs for the right hand side
   a1_[idt]->FormLinearSystem(fields_->dbc_dofs_, *fields_->dedt_, *rhs_, *A1_[idt], dEdt, *RHS_);

   // Solve for the time derivative of the electric field (true DoFs)
   pcg_[idt]->Mult(*RHS_, dEdt);

   // Distribute shared DoFs to relevant processors
   a1_[idt]->RecoverFEMSolution(dEdt, *rhs_, *fields_->dedt_);
}

double
VMSolver::GetMaximumTimeStep() const
{
   if ( dtMax_ > 0.0 )
   {
      return dtMax_;
   }

   mfem::HypreParVector * v0 = new mfem::HypreParVector(fields_->HCurlFESpace_);
   mfem::HypreParVector * v1 = new mfem::HypreParVector(fields_->HCurlFESpace_);
   mfem::HypreParVector * u0 = new mfem::HypreParVector(fields_->HDivFESpace_);

   v0->Randomize(1234);

   int iter = 0, nstep = 20;
   double dt0 = 1.0, dt1 = 1.0, change = 1.0, ptol = 0.001;

   // Create Solver assuming no loss operators
   setupSolver(0, 0.0);

   // Use power method to approximate the largest eigenvalue of the update
   // operator.
   while ( iter < nstep && change > ptol )
   {
      double normV0 = InnerProduct(*v0,*v0);
      *v0 /= sqrt(normV0);

      NegCurl_->Mult(*v0,*u0);
      M2MuInv_->Mult(*u0,*HD_);
      NegCurl_->MultTranspose(*HD_,*RHS_);

      pcg_[0]->Mult(*RHS_,*v1);

      double lambda = InnerProduct(*v0,*v1);
      dt1 = 2.0/sqrt(lambda);
      change = fabs((dt1-dt0)/dt0);
      dt0 = dt1;

      if ( fields_->myid_ == 0 && fields_->logging_ > 1 )
      {
         std::cout << iter << ":  " << dt0 << " " << change << std::endl;
      }

      std::swap(v0, v1);

      iter++;
   }

   delete v0;
   delete v1;
   delete u0;

   return dt0;
}

double
VMSolver::GetEnergy() const
{
   double energy = 0.0;

   A1_[0]->Mult(*fields_->E_,*RHS_);
   M2MuInv_->Mult(*fields_->B_,*HD_);

   energy = param_->epsilon0*InnerProduct(*fields_->E_,*RHS_) + InnerProduct(*fields_->B_,*HD_)/(param_->N_curlB*param_->mu0);

   return 0.5 * energy;
}

#endif // MFEM_USE_MPI
