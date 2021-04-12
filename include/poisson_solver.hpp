
#ifndef POISSON_SOLVER
#define POISSON_SOLVER

#include "field_functions.hpp"
#include "pfem_extras.hpp"
#include "mesh_extras.hpp"

#ifdef MFEM_USE_MPI

#include <string>
#include <map>

class PoissonSolver
{
public:
   PoissonSolver(ParametersInput & param, VMFields & fields, VMParticlesArray & particles);
   ~PoissonSolver();

   HYPRE_Int GetProblemSize();

   void PrintSizes();

   void Assemble();

   void Update();

   void Solve();

   void GetErrorEstimates(mfem::Vector & errors);

   const mfem::ParGridFunction & GetVectorPotential() { return *phi_; }

private:

   ParametersInput * param_;
   VMFields * fields_;
   VMParticlesArray * particles_;

   mfem::Vector     * dbcv_; // Corresponding Dirichlet Values

   mfem::ParBilinearForm * divEpsGrad_; // Laplacian operator

   mfem::ParLinearForm * rhod_; // Dual of Volumetric Charge Density Source

   mfem::common::ParDiscreteGradOperator * grad_; // For Computing E from phi
   mfem::common::ParDiscreteDivOperator  * div_;  // For Computing rho from D

   mfem::ParGridFunction * phi_;       // Electric Scalar Potential

   mfem::ConstantCoefficient oneCoef_;   // Coefficient equal to 1

   std::vector<mfem::DeltaCoefficient*> point_charges_;

   mfem::Array<int> ess_bdr_tdofs_; // Essential Boundary Condition DoFs
};

#endif // MFEM_USE_MPI

#endif // POISSON_SOLVER
