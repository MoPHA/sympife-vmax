
#include "field_functions.hpp"


#ifdef MFEM_USE_MPI


VMFields::VMFields(ParametersInput & param, mfem::ParMesh & pmesh)
   : myid_(0),
     num_procs_(1),
     logging_(1),
     param_(&param),
     pmesh_(&pmesh),
     face2elem_(pmesh.GetFaceToAllElementTable()),
     e_(NULL),
     b_(NULL),
     rho_(NULL),
     j_(NULL),
     dedt_(NULL),
     jd_(NULL),
     E_(NULL),
     B_(NULL),
     HCurlFESpace_(NULL),
     HDivFESpace_(NULL),
     epsCoef_(NULL),
     curlBCoef_(NULL),
     sigmaCoef_(NULL),
     etaInvCoef_(NULL),
     eCoef_(NULL),
     bCoef_(NULL),
     rhoCoef_(NULL),
     jCoef_(NULL),
     dEdtBCCoef_(NULL),
     visit_dc_(NULL),
     paraview_dc_(NULL)
{
   // Initialize MPI variables
   MPI_Comm_size(pmesh_->GetComm(), &num_procs_);
   MPI_Comm_rank(pmesh_->GetComm(), &myid_);


   // Define compatible parallel finite element spaces on the parallel
   // mesh. Here we use arbitrary order H1, Nedelec, and Raviart-Thomas finite
   // elements.
   H1FESpace_ = new mfem::common::H1_ParFESpace(pmesh_,param_->fe_order,pmesh_->Dimension());
   HCurlFESpace_ = new mfem::common::ND_ParFESpace(pmesh_,param_->fe_order,pmesh_->Dimension());
   HDivFESpace_  = new mfem::common::RT_ParFESpace(pmesh_,param_->fe_order,pmesh_->Dimension());
   L2FESpace_  = new mfem::common::L2_ParFESpace(pmesh_,param_->fe_order-1,pmesh_->Dimension());

   // Check for absorbing materials or boundaries
   lossy_ = param_->bc_absorbing.Size() > 0 ||  param_->sigma_sphere.Size() > 0 ;

   // Build grid functions
   e_    = new mfem::ParGridFunction(HCurlFESpace_);
   dedt_ = new mfem::ParGridFunction(HCurlFESpace_);
   b_    = new mfem::ParGridFunction(HDivFESpace_);
   rho_    = new mfem::ParGridFunction(L2FESpace_);
   *rho_ = 0.0;

   E_ = e_->ParallelProject();
   B_ = b_->ParallelProject();

   SetMaterialCoefficients();

   SetBoundaryConditions();

   SetInitialFields();

   // Initialize dedt to zero
   *dedt_ = 0.0;

}

HYPRE_Int
VMFields::GetProblemSize()
{
   return HCurlFESpace_->GlobalTrueVSize();
}

void
VMFields::PrintSizes()
{
   HYPRE_Int size_nd = HCurlFESpace_->GlobalTrueVSize();
   HYPRE_Int size_rt = HDivFESpace_->GlobalTrueVSize();
   if ( myid_ == 0 )
   {
      std::cout << "Number of H(Curl) unknowns: " << size_nd << std::endl;
      std::cout << "Number of H(Div)  unknowns: " << size_rt << std::endl << std::flush;
   }
}

void
VMFields::SetMaterialCoefficients()
{
   // Electric permittivity
   if ( param_->eps_sphere.Size() > 0 )
   {
      if ( myid_ == 0 && logging_ > 0 )
      {
         std::cout << "Creating Permittivity Coefficient" << std::endl;
      }
      epsCoef_ = new mfem::FunctionCoefficient([this](const mfem::Vector &x) { return epsilon(x); });
   }
   else
   {
      epsCoef_ = new mfem::ConstantCoefficient(1.0);
   }

   // Inverse of the magnetic permeability
   if ( param_->mu_shell.Size() > 0 )
   {
      if ( myid_ == 0 && logging_ > 0 )
      {
         std::cout << "Creating Permeability Coefficient" << std::endl;
      }
      mfem::FunctionCoefficient muInvCoef_([this](const mfem::Vector &x) { return muInv(x); });
      curlBCoef_ = new mfem::ProductCoefficient(param_->N_curlB, muInvCoef_);
   }
   else
   {
      curlBCoef_ = new mfem::ConstantCoefficient(param_->N_curlB);
   }

   // Electric conductivity
   if ( param_->sigma_sphere.Size() > 0 )
   {
      if ( myid_ == 0 && logging_ > 0 )
      {
         std::cout << "Creating Conductivity Coefficient" << std::endl;
      }
      sigmaCoef_ = new mfem::FunctionCoefficient([this](const mfem::Vector &x) { return sigma(x); });
   }

   // Current source
   if ( param_->dipole_pulse.Size() > 0 )
   {
      if ( myid_ == 0 && logging_ > 0 )
      {
         std::cout << "Creating Current Source" << std::endl;
      }

      jCoef_ = new mfem::VectorFunctionCoefficient(3,
         [this](const mfem::Vector &x, double t, mfem::Vector &j)
         { return j_src(x, t, j); });

      j_  = new mfem::ParGridFunction(HCurlFESpace_);
      j_->ProjectCoefficient(*jCoef_);

      jd_ = new mfem::ParLinearForm(HCurlFESpace_);
      jd_->AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*jCoef_));
      jd_->Assemble();
   }

   // Point particles
   if ( param_->point_charges.Size() > 0 )
   {
      if ( myid_ == 0 && logging_ > 0 )
      {
         std::cout << "Creating charge distribution" << std::endl;
      }
      int dim = pmesh_->Dimension();
      int npts = param_->point_charges.Size() / (dim + 1);

      mfem::DenseMatrix 	point_mat(dim,npts);
      for (int i=0; i<npts; i++)
      {
         mfem::Vector cent(dim);
         for (int d=0; d<dim; d++)
         {
            cent[d] = param_->point_charges[(dim + 1) * i + d];
         }
         point_mat.SetCol(i,cent);
      }

      mfem::Array<int>  	elem_ids;
      mfem::Array<mfem::IntegrationPoint>  	ips;
      pmesh_->FindPoints(point_mat, elem_ids, ips);

      for (int i=0; i<npts; i++)
      {
         if ( elem_ids[i] > 0 )
         {
            double vol = pmesh_->GetElementVolume(elem_ids[i]);

            mfem::Vector vals;
            mfem::Array<int> dofs;

            const mfem::FiniteElement *fe = rho_->FESpace()->GetFE(elem_ids[i]);
            vals.SetSize(fe->GetDof());
            dofs.SetSize(fe->GetDof());
            fe->CalcShape(ips[i], vals);
            rho_->FESpace()->GetElementDofs(elem_ids[i], dofs);
            vals *= param_->point_charges[(dim + 1) * i + dim]/vol;
            rho_->AddElementVector(dofs,vals);
         }
      }
   }
}




void
VMFields::SetBoundaryConditions()
{
   // Impedance of free space
   if ( param_->bc_absorbing.Size() > 0 )
   {
      if ( myid_ == 0 && logging_ > 0 )
      {
         std::cout << "Creating Admittance Coefficient" << std::endl;
      }

      mfem::common::AttrToMarker(pmesh_->bdr_attributes.Max(), param_->bc_absorbing, abc_marker_);
      etaInvCoef_ = new mfem::ConstantCoefficient(1.0/param_->N_eta);
   }

   // Electric Field Boundary Condition
   if ( param_->bc_dirichlet.Size() > 0 )
   {
      if ( myid_ == 0 && logging_ > 0 )
      {
         std::cout << "Configuring Dirichlet BC" << std::endl;
      }

      mfem::common::AttrToMarker(pmesh_->bdr_attributes.Max(), param_->bc_dirichlet, dbc_marker_);
      HCurlFESpace_->GetEssentialTrueDofs(dbc_marker_, dbc_dofs_);

      dEdtBCCoef_ = new mfem::VectorFunctionCoefficient(3,
            [this](const mfem::Vector &x, double t, mfem::Vector &dE) { return dEdt_bc(x, t, dE); });

   }
}

void
VMFields::SetInitialFields()
{

   eCoef_ = new mfem::VectorFunctionCoefficient(3,
      [this](const mfem::Vector &x, mfem::Vector &E) { return E_init(x, E); });

   e_->ProjectCoefficient(*eCoef_);
   e_->ParallelProject(*E_);

   bCoef_ = new mfem::VectorFunctionCoefficient(3,
      [this](const mfem::Vector &x, mfem::Vector &B) { return B_init(x, B); });

   b_->ProjectCoefficient(*bCoef_);
   b_->ParallelProject(*B_);
}

void
VMFields::E_function(const mfem::Vector &x, mfem::Vector &E)
{
   E.SetSize(3);
   E = 0.0;
}

void
VMFields::B_function(const mfem::Vector &x, mfem::Vector &B)
{
   B.SetSize(3);
   B = 0.0;
}

void
VMFields::SyncGridFuncs()
{
   e_->Distribute(*E_);
   b_->Distribute(*B_);
}


// A sphere with constant permittivity.  The sphere has a radius, center, and
// permittivity specified on the command line and stored in ds_params_.
double
VMFields::dielectric_sphere(const mfem::Vector &x)
{
   double r2 = 0.0;

   for (int i=0; i<x.Size(); i++)
   {
      r2 += (x(i)-param_->eps_sphere[i])*(x(i)-param_->eps_sphere[i]);
   }

   if ( sqrt(r2) <= param_->eps_sphere[x.Size()] )
   {
      return param_->eps_sphere[x.Size()+1];
   }
   return 1.0;
}

// A spherical shell with constant permeability.  The sphere has inner and outer
// radii, center, and relative permeability specified on the command line and
// stored in ms_params_.
double
VMFields::magnetic_shell(const mfem::Vector &x)
{
   double r2 = 0.0;

   for (int i=0; i<x.Size(); i++)
   {
      r2 += (x(i)-param_->mu_shell[i])*(x(i)-param_->mu_shell[i]);
   }

   if ( sqrt(r2) >= param_->mu_shell[x.Size()] &&
        sqrt(r2) <= param_->mu_shell[x.Size()+1] )
   {
      return param_->mu_shell[x.Size()+2];
   }
   return 1.0;
}

// A sphere with constant conductivity.  The sphere has a radius, center, and
// conductivity specified on the command line and stored in ls_params_.
double
VMFields::conductive_sphere(const mfem::Vector &x)
{
   double r2 = 0.0;

   for (int i=0; i<x.Size(); i++)
   {
      r2 += (x(i)-param_->sigma_sphere[i])*(x(i)-param_->sigma_sphere[i]);
   }

   if ( sqrt(r2) >= param_->sigma_sphere[x.Size()] )
   {
      return param_->sigma_sphere[x.Size()+1];
   }
   return 0.0;
}

// A cylindrical rod of current density.  The rod has two axis end points, a
// radius, a current amplitude in Amperes, a center time, and a width.  All of
// these parameters are stored in dipole_pulse.
void
VMFields::dipole_function(const mfem::Vector &x, double t, mfem::Vector &j)
{
   MFEM_ASSERT(x.Size() == 3, "current source requires 3D space.");

   j.SetSize(x.Size());
   j = 0.0;

   mfem::Vector  v(x.Size());  // Normalized Axis vector
   mfem::Vector xu(x.Size());  // x vector relative to the axis end-point

   xu = x;

   for (int i=0; i<x.Size(); i++)
   {
      xu[i] -= param_->dipole_pulse[i];
      v[i]   = param_->dipole_pulse[x.Size()+i] - param_->dipole_pulse[i];
   }

   double h = v.Norml2();

   if ( h == 0.0 )
   {
      return;
   }
   v /= h;

   double r = param_->dipole_pulse[2*x.Size()+0];
   double a = param_->dipole_pulse[2*x.Size()+1];
   double b = param_->dipole_pulse[2*x.Size()+2];
   double c = param_->dipole_pulse[2*x.Size()+3];

   double xv = xu * v;

   // Compute perpendicular vector from axis to x
   xu.Add(-xv, v);

   double xp = xu.Norml2();

   if ( xv >= 0.0 && xv <= h && xp <= r )
   {
      j = v;
   }

   j *= a * (t - b) * std::exp(-0.5 * std::pow((t-b)/c, 2)) / (c * c);
}

void
VMFields::dEdtBC_function(const mfem::Vector &x, double t, mfem::Vector &dE)
{
   dE.SetSize(3);
   dE = 0.0;
}

void
VMFields::RegisterVisItFields(mfem::VisItDataCollection & visit_dc)
{
   visit_dc_ = &visit_dc;

   visit_dc.RegisterField("E", e_);
   visit_dc.RegisterField("B", b_);
   visit_dc.RegisterField("rho", rho_);
   if ( j_ )
   {
      visit_dc.RegisterField("J", j_);
   }
}

void
VMFields::WriteVisItFields(int it, double t)
{
   if ( visit_dc_ )
   {
      if ( myid_ == 0 && logging_ > 1 )
      { std::cout << "Writing VisIt files ..." << std::flush; }

      if ( j_ )
      {
         jCoef_->SetTime(t);
         j_->ProjectCoefficient(*jCoef_);
      }

      visit_dc_->SetCycle(it);
      visit_dc_->SetTime(t);
      visit_dc_->Save();

      if ( myid_ == 0 && logging_ > 1 ) { std::cout << " " << std::endl << std::flush; }
   }
}

void
VMFields::RegisterParaViewFields(mfem::ParaViewDataCollection & paraview_dc)
{
   paraview_dc_ = &paraview_dc;

   paraview_dc.SetPrefixPath("ParaView");

   paraview_dc.RegisterField("E", e_);
   paraview_dc.RegisterField("B", b_);
   paraview_dc.RegisterField("rho", rho_);
   if ( j_ )
   {
      paraview_dc.RegisterField("J", j_);
   }
}

void
VMFields::WriteParaViewFields(int it, double t)
{
   if ( paraview_dc_ )
   {
      if ( myid_ == 0 && logging_ > 1 )
      { std::cout << "Writing ParaView files ..." << std::flush; }

      paraview_dc_->SetLevelsOfDetail(param_->fe_order);
      paraview_dc_->SetDataFormat(mfem::VTKFormat::BINARY);
      paraview_dc_->SetHighOrderOutput(true);

      if ( j_ )
      {
         jCoef_->SetTime(t);
         j_->ProjectCoefficient(*jCoef_);
      }

      paraview_dc_->SetCycle(it);
      paraview_dc_->SetTime(t);
      paraview_dc_->Save();

      if ( myid_ == 0 && logging_ > 1 ) { std::cout << " " << std::endl << std::flush; }
   }
}

#endif // MFEM_USE_MPI
