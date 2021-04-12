
#ifndef FIELD_FUNCTIONS
#define FIELD_FUNCTIONS


#include "parameters.hpp"
#include "particle_functions.hpp"
#include "pfem_extras.hpp"
#include "mesh_extras.hpp"


#ifdef MFEM_USE_MPI


#include <string>
#include <map>


class VMFields
{
public:
   VMFields(ParametersInput & param, mfem::ParMesh & pmesh);

   int GetLogging() const { return logging_; }
   void SetLogging(int logging) { logging_ = logging; }

   double dielectric_sphere(const mfem::Vector &);
   double epsilon(const mfem::Vector &x) { return dielectric_sphere(x); }

   double magnetic_shell(const mfem::Vector &);
   double muInv(const mfem::Vector & x) { return 1.0/magnetic_shell(x); }

   double conductive_sphere(const mfem::Vector &);
   double sigma(const mfem::Vector &x) { return conductive_sphere(x); }

   void dipole_function(const mfem::Vector &x, double t, mfem::Vector &j);
   void j_src(const mfem::Vector &x, double t, mfem::Vector &j) { dipole_function(x, t, j); }

   void dEdtBC_function(const mfem::Vector &x, double t, mfem::Vector &dE);
   void dEdt_bc(const mfem::Vector &x, double t, mfem::Vector &dE) { dEdtBC_function(x, t, dE); }

   void E_function(const mfem::Vector &x, mfem::Vector &E);
   void E_init(const mfem::Vector &x, mfem::Vector &E) { E_function(x, E); }

   void B_function(const mfem::Vector &x, mfem::Vector &B);
   void B_init(const mfem::Vector &x, mfem::Vector &B) { B_function(x, B); }

   HYPRE_Int GetProblemSize();

   void PrintSizes();

   void SetMaterialCoefficients();
   void SetBoundaryConditions();

   void SetInitialFields();

   mfem::Vector & GetEField() { return *E_; }
   mfem::Vector & GetBField() { return *B_; }

   void SyncGridFuncs();

   void RegisterVisItFields(mfem::VisItDataCollection & visit_dc);

   void WriteVisItFields(int it = 0, double t = 0);

   void RegisterParaViewFields(mfem::ParaViewDataCollection & paraview_dc);

   void WriteParaViewFields(int it = 0, double t = 0);


   int myid_;
   int num_procs_;
   int logging_;

   bool lossy_;

   ParametersInput * param_;

   mfem::ParMesh * pmesh_;

   mfem::Table * face2elem_;

   mfem::common::H1_ParFESpace * H1FESpace_;    // Continuous space for phi
   mfem::common::ND_ParFESpace * HCurlFESpace_; // Tangentially continuous space for E
   mfem::common::RT_ParFESpace * HDivFESpace_;  // Normally continuous space for D
   mfem::common::L2_ParFESpace * L2FESpace_;    // Discontinuous space for rho

   mfem::ParGridFunction * e_;    // Electric Field (HCurl)
   mfem::ParGridFunction * b_;    // Magnetic Flux (HDiv)
   mfem::ParGridFunction * rho_;  // Charge density (L2)
   mfem::ParGridFunction * j_;    // Volumetric Current Density (HCurl)
   mfem::ParGridFunction * dedt_; // Time Derivative of Electric Field (HCurl)
   mfem::ParLinearForm   * jd_;   // Dual of current density (HCurl)

   mfem::HypreParVector * E_; // Current value of the electric field DoFs
   mfem::HypreParVector * B_; // Current value of the magnetic flux DoFs

   mfem::Coefficient       * epsCoef_;    // Electric Permittivity Coefficient
   mfem::Coefficient       * curlBCoef_;  // Magnetic Permeability Coefficient
   mfem::Coefficient       * sigmaCoef_;  // Electric Conductivity Coefficient
   mfem::Coefficient       * etaInvCoef_; // Admittance Coefficient
   mfem::DeltaCoefficient  * rhoCoef_;    // Initial charge density
   mfem::VectorCoefficient * eCoef_;      // Initial Electric Field
   mfem::VectorCoefficient * bCoef_;      // Initial Magnetic Flux
   mfem::VectorCoefficient * jCoef_;      // Time dependent current density
   mfem::VectorCoefficient * dEdtBCCoef_; // Time dependent boundary condition

   double (*eps_    )(const mfem::Vector&);
   double (*muInv_  )(const mfem::Vector&);

   // Array of 0's and 1's marking the location of absorbing surfaces
   mfem::Array<int> abc_marker_;

   // Array of 0's and 1's marking the location of Dirichlet boundaries
   mfem::Array<int> dbc_marker_;

   // Dirichlet degrees of freedom
   mfem::Array<int>   dbc_dofs_;

   std::vector<mfem::DeltaCoefficient*> point_charges_;

   private:

   // Data collection used to write VisIt files
   mfem::VisItDataCollection * visit_dc_;

   // Data collection used to write ParaView files
   mfem::ParaViewDataCollection * paraview_dc_;
};



#endif // MFEM_USE_MPI

#endif // FIELD_FUNCTIONS
