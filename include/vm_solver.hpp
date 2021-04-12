
#ifndef VM_SOLVER
#define VM_SOLVER


#include "pfem_extras.hpp"
#include "mesh_extras.hpp"
#include "field_functions.hpp"


#ifdef MFEM_USE_MPI


#include <string>
#include <map>


class VMSolver : public mfem::TimeDependentOperator
{
public:
   VMSolver(ParametersInput & param, VMFields & fields, VMParticlesArray & particles);

   ~VMSolver();

   HYPRE_Int GetProblemSize();

   void PrintSizes();

   void Mult(const mfem::Vector &B, mfem::Vector &dEdt) const;

   void ImplicitSolve(const double dt, const mfem::Vector &x, mfem::Vector &k);

   double GetMaximumTimeStep() const;

   double GetEnergy() const;

   mfem::Operator & GetNegCurl() { return *NegCurl_; }

   ParametersInput * param_;
   VMFields * fields_;
   VMParticlesArray * particles_;

private:

   // This method alters mutable member data
   void setupSolver(const int idt, const double dt) const;

   void implicitSolve(const double dt, const mfem::Vector &x, mfem::Vector &k) const;


   bool lossy_;

   double dtMax_;   // Maximum stable time step
   double dtScale_; // Used to scale dt before converting to an integer

   mfem::ParBilinearForm * hDivMassMuInv_;
   mfem::ParBilinearForm * hCurlLosses_;
   mfem::ParMixedBilinearForm * weakCurlMuInv_;

   mfem::common::ParDiscreteCurlOperator * Curl_;

   mfem::ParGridFunction * rhs_;  // Dual of displacement current, rhs vector (HCurl)

   mfem::HypreParMatrix * M1Losses_;
   mfem::HypreParMatrix * M2MuInv_;
   mfem::HypreParMatrix * NegCurl_;
   mfem::HypreParMatrix * WeakCurlMuInv_;
   mutable mfem::HypreParVector * HD_; // Used in energy calculation
   mutable mfem::HypreParVector * RHS_;

   // High order symplectic integration requires partial time steps of differing
   // lengths. If losses are present the system matrix includes a portion scaled
   // by the time step. Consequently, high order time integration requires
   // different system matrices. The following maps contain various objects that
   // depend on the time step.
   mutable std::map<int, mfem::ParBilinearForm *> a1_;
   mutable std::map<int, mfem::HypreParMatrix  *> A1_;
   mutable std::map<int, mfem::Coefficient     *> dtCoef_;
   mutable std::map<int, mfem::Coefficient     *> dtSigmaCoef_;
   mutable std::map<int, mfem::Coefficient     *> dtEtaInvCoef_;
   mutable std::map<int, mfem::HypreDiagScale  *> diagScale_;
   mutable std::map<int, mfem::HyprePCG        *> pcg_;
};



#endif // MFEM_USE_MPI

#endif // MFEM_MAXWELL_SOLVER
