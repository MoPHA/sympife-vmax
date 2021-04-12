
#include "particle_functions.hpp"


#ifdef MFEM_USE_MPI

static std::random_device rd;
static std::mt19937 rng(rd());


VMParticlesArray::VMParticlesArray(ParametersInput & param, mfem::ParMesh & pmesh)
   : myid_(0),
   num_procs_(0),
   param_(&param),
   pmesh_(&pmesh)

{
   // Initialize MPI variables
   MPI_Comm_size(pmesh_->GetComm(), &num_procs_);
   MPI_Comm_rank(pmesh_->GetComm(), &myid_);
   std::srand((unsigned) std::time(0)*myid_);

   particlesArray_ = new mfem::Array2D<VMParticles*>(pmesh_->GetNE(),param_->num_species);

   for (int i=0; i<pmesh_->GetNE(); i++)
   {
      for (int j=0; j<param_->num_species; j++)
      {

         int number = 1000;
         int ij = param_->num_species*i+j;
         (*particlesArray_)(i,j) = new VMParticles(number,1.0/sqrt(2.0));

      }

   }

}


double
VMParticlesArray::n_init(const mfem::Vector &)
{

   //return (cos(2.0*M_PI*x[0])*cos(2.0*M_PI*x[1])*cos(2.0*M_PI*x[2]))+2.0;
   return param_->rho0;

}


VMParticles::VMParticles(int number,double velocity)
   : size(number)
{


   x1_ = new mfem::Vector(number);
   x1_->Randomize(std::rand()%10000);
   x2_ = new mfem::Vector(number);
   x2_->Randomize(std::rand()%10000);
   x3_ = new mfem::Vector(number);
   x3_->Randomize(std::rand()%10000);



   v1_ = new mfem::Vector(number);
   v2_ = new mfem::Vector(number);
   v3_ = new mfem::Vector(number);

   std::default_random_engine generator;
   std::normal_distribution<double> distribution(0.0,velocity);

   for (int i=0; i<number; i++)
   {
      (*v1_)(i) = distribution(rng);
      (*v2_)(i) = distribution(rng);
      (*v3_)(i) = distribution(rng);
   }
}

#endif // MFEM_USE_MPI
