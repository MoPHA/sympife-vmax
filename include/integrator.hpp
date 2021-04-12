
#ifndef SYM_INTEGRATOR
#define SYM_INTEGRATOR

#include "pfem_extras.hpp"
#include "vm_solver.hpp"

#ifdef MFEM_USE_MPI


#include <string>
#include <map>



// Explicit Hamiltonian Splitting PIC symplectic algorithm
class EHSMaxwellSolver
{
public:
   EHSMaxwellSolver() : F_(NULL), P_(NULL) {}

   virtual void Init(mfem::Operator &P, mfem::TimeDependentOperator & F);

   virtual void Step(mfem::Vector &B, mfem::Vector &E, double &t, double &dt) = 0;

   virtual void Run(mfem::Vector &B, mfem::Vector &E, double &t, double &dt, double tf)
   {
      //initial full step with first alpha/2
      //steps with first alpha
      while (t < tf) { Step(B, E, t, dt); }
      //final thetax with alpha/2
   }

   virtual ~EHSMaxwellSolver() {}

protected:
   mfem::TimeDependentOperator * F_; // E_{i+1} = E_{i} + dt F(B_{i})
   mfem::Operator              * P_; // B_{i+1} = B_{i} + dt P(E_{i+1})

   mutable mfem::Vector dE_;
   mutable mfem::Vector dB_;
};

/// Variable order Symplectic Integration Algorithm (orders 1-4)
class EHSVMaxwellSolver : public EHSMaxwellSolver
{
public:
   EHSVMaxwellSolver(int order);
   void Step(mfem::Vector &B, mfem::Vector &E, double &t, double &dt) override;

private:
   int order_;

   mfem::Array<double> a_;
   mfem::Array<double> b_;
};





/// Variable order Symplectic Integration Algorithm (orders 1-4)
class EHSVSolver
{
public:
   EHSVSolver(int order, VMSolver &vmsolver, mfem::Operator &P, mfem::TimeDependentOperator & F);

   void Step1(mfem::Vector &B, mfem::Vector &E, double &t, double &dt);

   void Step2(mfem::Vector &B, mfem::Vector &E, double &t, double &dt);

   void Step4(mfem::Vector &B, mfem::Vector &E, double &t, double &dt);

   void Run(mfem::Vector &B, mfem::Vector &E, double &t, double &dt, double tf);


   void Push_x1(double dt);
   void Push_x2(double dt);
   void Push_x3(double dt);

   ~EHSVSolver() {}

   VMSolver * vmsolver_;

private:
   int order_, blocks_;
   double alpha_, beta_;
   mfem::Array<double> a_;

protected:
   mfem::TimeDependentOperator * F_; // E_{i+1} = E_{i} + dt F(B_{i})
   mfem::Operator              * P_; // B_{i+1} = B_{i} + dt P(E_{i+1})

   mutable mfem::Vector dE_;
   mutable mfem::Vector dB_;
};



#endif //SYM_INTEGRATOR

#endif // MFEM_USE_MPI
