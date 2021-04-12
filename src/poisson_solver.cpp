
#include "poisson_solver.hpp"

#ifdef MFEM_USE_MPI

PoissonSolver::PoissonSolver(ParametersInput & param, VMFields & fields, VMParticlesArray & particles)
   : fields_(&fields),
     param_(&param),
     particles_(&particles),
     divEpsGrad_(NULL),
     rhod_(NULL),
     grad_(NULL),
     phi_(NULL),
     oneCoef_(1.0),
     point_charges_(0)
{

   // Bilinear Forms
   divEpsGrad_  = new mfem::ParBilinearForm(fields_->H1FESpace_);
   divEpsGrad_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*(fields_->epsCoef_)));


   rhod_   = new mfem::ParLinearForm(fields_->H1FESpace_);

   // Discrete derivative operator
   grad_ = new mfem::common::ParDiscreteGradOperator(fields_->H1FESpace_, fields_->HCurlFESpace_);

   // Build grid functions
   phi_  = new mfem::ParGridFunction(fields_->H1FESpace_);


   if ( param_->point_charges.Size() > 0 )
   {
      int dim = fields_->pmesh_->Dimension();
      int npts = param_->point_charges.Size() / (dim + 1);
      point_charges_.resize(npts);

      mfem::Vector cent(dim);
      for (int i=0; i<npts; i++)
      {
         for (int d=0; d<dim; d++)
         {
            cent[d] = param_->point_charges[(dim + 1) * i + d];
         }
         double s = param_->point_charges[(dim + 1) * i + dim];

         point_charges_[i] = new mfem::DeltaCoefficient();
         point_charges_[i]->SetScale(s);
         point_charges_[i]->SetDeltaCenter(cent);

         rhod_->AddDomainIntegrator(new mfem::DomainLFIntegrator(*point_charges_[i]));
      }
   }
   else
   {
      int dim = fields_->pmesh_->Dimension();

      mfem::DeltaCoefficient* point_charge_;
      for (int i=0; i<fields_->pmesh_->GetNE(); i++)
      {
         for (int j=0; j<param_->num_species; j++)
         {

            VMParticles* partsInElem = (*particles_->particlesArray_)(i,j);

            int npts = partsInElem->size;

            mfem::ElementTransformation * etrans = fields_->H1FESpace_->GetElementTransformation(i);

            mfem::Vector cent(dim);
            for (int k=0; k<npts; k++)
            {
               mfem::IntegrationPoint ip;
               ip.Set3((*partsInElem->x1_)(k),(*partsInElem->x2_)(k),(*partsInElem->x3_)(k));
               etrans->Transform(ip, cent);

               point_charge_ = new mfem::DeltaCoefficient();

               point_charge_->SetScale(param_->N_q_s[j]);
               //Doesn't seem to take the charge sign
               point_charge_->SetDeltaCenter(cent);

               rhod_->AddDomainIntegrator(new mfem::DomainLFIntegrator(*point_charge_));
            }
         }
      }
   }
}

PoissonSolver::~PoissonSolver()
{

   delete phi_;
   delete rhod_;

   delete grad_;

   delete divEpsGrad_;

   for (unsigned int i=0; i<point_charges_.size(); i++)
   {
      delete point_charges_[i];
   }
}

HYPRE_Int
PoissonSolver::GetProblemSize()
{
   return fields_->H1FESpace_->GlobalTrueVSize();
}

void
PoissonSolver::PrintSizes()
{
   HYPRE_Int size_h1 = fields_->H1FESpace_->GlobalTrueVSize();
   HYPRE_Int size_nd = fields_->HCurlFESpace_->GlobalTrueVSize();
   HYPRE_Int size_rt = fields_->HDivFESpace_->GlobalTrueVSize();
   HYPRE_Int size_l2 = fields_->L2FESpace_->GlobalTrueVSize();
   if (fields_->myid_ == 0)
   {
      std::cout << "Number of H1      unknowns: " << size_h1 << std::endl;
      std::cout << "Number of H(Curl) unknowns: " << size_nd << std::endl;
      std::cout << "Number of H(Div)  unknowns: " << size_rt << std::endl;
      std::cout << "Number of L2      unknowns: " << size_l2 << std::endl;
   }
}

void PoissonSolver::Assemble()
{
   if (fields_->myid_ == 0) { std::cout << "Assembling ... " << std::flush; }

   divEpsGrad_->Assemble();
   divEpsGrad_->Finalize();

   *rhod_ = 0.0;
   rhod_->Assemble();

   grad_->Assemble();
   grad_->Finalize();

   if (fields_->myid_ == 0) { std::cout << "done." << std::endl << std::flush; }
}

void
PoissonSolver::Update()
{
   if (fields_->myid_ == 0) { std::cout << "Updating ..." << std::endl; }

   // Inform the spaces that the mesh has changed
   // Note: we don't need to interpolate any GridFunctions on the new mesh
   // so we pass 'false' to skip creation of any transformation matrices.
   fields_->H1FESpace_->Update(false);
   fields_->HCurlFESpace_->Update(false);
   fields_->HDivFESpace_->Update(false);
   fields_->L2FESpace_->Update(false);

   // Inform the grid functions that the space has changed.
   phi_->Update();
   rhod_->Update();
   fields_->e_->Update();

   // Inform the bilinear forms that the space has changed.
   divEpsGrad_->Update();

   // Inform the other objects that the space has changed.
   grad_->Update();
}

void
PoissonSolver::Solve()
{
   if (fields_->myid_ == 0) { std::cout << "Running solver ... " << std::endl; }

   // Initialize the electric potential with its boundary conditions
   *phi_ = 0.0;

   if ( param_->bc_dirichlet.Size() > 0 )
   {
      // Apply piecewise constant boundary condition
      mfem::Array<int> dbc_bdr_attr(fields_->pmesh_->bdr_attributes.Max());
      for (int i=0; i<param_->bc_dirichlet.Size(); i++)
      {
         mfem::ConstantCoefficient voltage(0.0);
         dbc_bdr_attr = 0;
         if (param_->bc_dirichlet[i] <= dbc_bdr_attr.Size())
         {
            dbc_bdr_attr[param_->bc_dirichlet[i]-1] = 1;
         }
         phi_->ProjectBdrCoefficient(voltage, dbc_bdr_attr);
      }
   }

   // Determine the essential BC degrees of freedom
   if ( param_->bc_dirichlet.Size() > 0 )
   {
      // From user supplied boundary attributes
      fields_->H1FESpace_->GetEssentialTrueDofs(fields_->dbc_marker_, ess_bdr_tdofs_);
   }
   else
   {
      // Use the first DoF on processor zero by default
      if ( fields_->myid_ == 0 )
      {
         ess_bdr_tdofs_.SetSize(1);
         ess_bdr_tdofs_[0] = 0;
      }
   }

   // Apply essential BC and form linear system
   mfem::HypreParMatrix DivEpsGrad;
   mfem::HypreParVector Phi(fields_->H1FESpace_);
   mfem::HypreParVector RHS(fields_->H1FESpace_);

   divEpsGrad_->FormLinearSystem(ess_bdr_tdofs_, *phi_, *rhod_, DivEpsGrad,
                                 Phi, RHS);

   // Define and apply a parallel PCG solver for AX=B with the AMG
   // preconditioner from hypre.
   mfem::HypreBoomerAMG amg(DivEpsGrad);
   mfem::HyprePCG pcg(DivEpsGrad);
   pcg.SetTol(1e-12);
   pcg.SetMaxIter(500);
   pcg.SetPrintLevel(0);
   pcg.SetPreconditioner(amg);
   pcg.Mult(RHS, Phi);

   // Extract the parallel grid function corresponding to the finite
   // element approximation Phi. This is the local solution on each
   // processor.
   divEpsGrad_->RecoverFEMSolution(Phi, *rhod_, *phi_);

   // Compute the negative Gradient of the solution vector.  This is
   // the magnetic field corresponding to the scalar potential
   // represented by phi.
   grad_->Mult(*phi_, *(fields_->e_)); *(fields_->e_) *= -1.0;

   if (fields_->myid_ == 0) { std::cout << "Poisson solver done. " << std::endl; }
}

void
PoissonSolver::GetErrorEstimates(mfem::Vector & errors)
{
   if (fields_->myid_ == 0) { std::cout << "Estimating Error ... " << std::flush; }

   // Space for the discontinuous (original) flux
   mfem::DiffusionIntegrator flux_integrator(*(fields_->epsCoef_));
   mfem::L2_FECollection flux_fec(param_->fe_order, fields_->pmesh_->Dimension());
   // ND_FECollection flux_fec(param_->fe_order, fields_->pmesh_->Dimension());
   mfem::ParFiniteElementSpace flux_fes(fields_->pmesh_, &flux_fec, fields_->pmesh_->SpaceDimension());

   // Space for the smoothed (conforming) flux
   double norm_p = 1;
   mfem::RT_FECollection smooth_flux_fec(param_->fe_order-1, fields_->pmesh_->Dimension());
   mfem::ParFiniteElementSpace smooth_flux_fes(fields_->pmesh_, &smooth_flux_fec);

   L2ZZErrorEstimator(flux_integrator, *phi_,
                      smooth_flux_fes, flux_fes, errors, norm_p);

   if (fields_->myid_ == 0) { std::cout << "done." << std::endl; }
}

#endif // MFEM_USE_MPI
