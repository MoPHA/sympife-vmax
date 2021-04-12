
#ifndef PARTICLE_FUNCTIONS
#define PARTICLE_FUNCTIONS


#include "parameters.hpp"


#ifdef MFEM_USE_MPI


#include <string>
#include <map>
#include <random>




class VMParticles
{
public:
   VMParticles(int number,double velocity);
   int size;

   mfem::Vector * x1_;
   mfem::Vector * x2_;
   mfem::Vector * x3_;
   mfem::Vector * v1_;
   mfem::Vector * v2_;
   mfem::Vector * v3_;
};



class VMParticlesArray
{
public:
   VMParticlesArray(ParametersInput & param_, mfem::ParMesh & pmesh);

   double n_init(const mfem::Vector &);
   double nCoef(const mfem::Vector &x) { return n_init(x); }

   int myid_;
   int num_procs_;


   ParametersInput * param_;

   mfem::ParMesh * pmesh_;

   mfem::Array2D<VMParticles*> * particlesArray_;
};

#endif // MFEM_USE_MPI

#endif // PARTICLE_FUNCTIONS
