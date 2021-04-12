// This miniapp solves a simple 6D full-wave Vlasov-Maxwell problem using the
// coupled, first-order equations:
//
//                 epsilon dE/dt = Curl 1/mu B - sigma E - J
//                         dB/dt = - Curl E
//
// The permittivity function is that of the vacuum with an optional dielectric
// sphere. The permeability function is that of the vacuum with an optional
// diamagnetic or paramagnetic spherical shell. The optional conductivity
// function is also a user-defined sphere.
//
// The optional current density is a pulse of current in the shape of a cylinder
// with a time dependence resembling the derivative of a Gaussian distribution.
//
// Boundary conditions can be 'natural' meaning zero tangential current,
// 'Dirichlet' which sets the time-derivative of the tangential components of E,
// or 'absorbing' (we use a simple Sommerfeld first order absorbing boundary
// condition).
//
// We discretize the electric field with H(Curl) finite elements (Nedelec edge
// elements) and the magnetic flux with H(Div) finite elements (Raviart-Thomas
// elements).
//
// The symplectic time integration algorithm used below is designed to conserve
// energy unless lossy materials or absorbing boundary conditions are used.
// When losses are expected, the algorithm uses an implicit method which
// includes the loss operators in the left hand side of the linear system.
//
// For increased accuracy the time integration order can be set to 2 or 4
// (the default is 1st order).
//
// Compile with: mkdir build ; cd build ; cmake .. ; cmake --build .
//
// Run with: mpirun -np 4 SymPiFE -c "relative/path/to/config_file.txt"
//
//   By default the config file path is "../config.txt"
//     mpirun -np 4 SymPiFE

#include "parameters.hpp"
#include "field_functions.hpp"
#include "particle_functions.hpp"
#include "poisson_solver.hpp"
#include "vm_solver.hpp"
#include "integrator.hpp"
#include <fstream>
#include <iostream>


// Scale factor between input time units and seconds
static double tScale_ = 1e-9;  // Input time in nanosecond

int SnapTimeStep(double tmax, double dtmax, double & dt);

int main(int argc, char *argv[])
{
   mfem::MPI_Session mpi(argc, argv);


   // Parse command-line options.
   const char *config_file = "../config.txt";

   //std::ofstream out("out.txt");
   //std::streambuf *coutbuf = std::cout.rdbuf();
   //std::cout.rdbuf(out.rdbuf());

   mfem::OptionsParser args(argc, argv);
   args.AddOption(&config_file, "-c", "--config",
   "Config file to use.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.Root())
      {
         args.PrintUsage(std::cout);
      }
      return 1;
   }
   if (mpi.Root())
   {
      args.PrintOptions(std::cout);
   }

   ParametersInput param(config_file);

   // Read the (serial) mesh from the given mesh file on all processors.  We can
   // handle triangular, quadrilateral, tetrahedral, hexahedral, surface and
   // volume meshes with the same code.
   mfem::Mesh *mesh;
   std::ifstream imesh(param.mesh_path);

   if (!imesh)
   {
      if (mpi.Root())
      {
         std::cerr << "\nCan not open mesh file: " << param.mesh_path << '\n' << std::endl;
      }
      return 2;
   }
   mesh = new mfem::Mesh(imesh, 1, 1);
   imesh.close();

   // Project a NURBS mesh to a piecewise-quadratic curved mesh
   if (mesh->NURBSext)
   {
      mesh->UniformRefinement();
      if (param.mesh_serial_ref > 0) { param.mesh_serial_ref--; }

      mesh->SetCurvature(2);
   }

   // Refine the serial mesh on all processors to increase the resolution. In
   // this example we do 'ref_levels' of uniform refinement.
   for (int l = 0; l < param.mesh_serial_ref; l++)
   {
      mesh->UniformRefinement();
   }

   // Define a parallel mesh by a partitioning of the serial mesh. Refine this
   // mesh further in parallel to increase the resolution. Once the parallel
   // mesh is defined, the serial mesh can be deleted.
   mfem::ParMesh pmesh(MPI_COMM_WORLD, *mesh, NULL);
   delete mesh;

   // Refine this mesh in parallel to increase the resolution.
   for (int l = 0; l < param.mesh_para_ref; l++)
   {
      pmesh.UniformRefinement();
   }

   // Create the Electromagnetic solver

   VMFields Fields(param, pmesh);

   if ( mpi.Root() )
   {
      std::cout << "VM fields created" << std::endl;
   }

   VMParticlesArray Particles(param, pmesh);

   /*PoissonSolver PoissonSolver(param, Fields, Particles);

   // Display the current number of DoFs in each finite element space
   PoissonSolver.PrintSizes();

   // Assemble all forms
   PoissonSolver.Assemble();

   // Solve the system and compute any auxiliary fields
   PoissonSolver.Solve();*/

   VMSolver VlasovMaxwell(param, Fields, Particles);

   if ( mpi.Root() )
   {
      std::cout << "VM solver created" << std::endl;
   }

   // Display the current number of DoFs in each finite element space
   Fields.PrintSizes();

   // Compute the energy of the initial fields
   double energy = VlasovMaxwell.GetEnergy();
   if ( mpi.Root() )
   {
      std::cout << "Energy(" << param.time_ini << "ns):  " << energy << "J" << std::endl;
   }

   // Approximate the largest stable time step
   double dtmax = VlasovMaxwell.GetMaximumTimeStep();

   // Convert times from nanoseconds to seconds
   param.time_ini *= tScale_;
   param.time_fin *= tScale_;
   param.time_diag *= tScale_;

   if ( mpi.Root() )
   {
      std::cout << "Maximum Time Step:     " << dtmax / tScale_ << "ns" << std::endl;
   }

   // Round down the time step so that tf-ti is an integer multiple of dt
   int nsteps = SnapTimeStep(param.time_fin-param.time_ini, param.dt_safe * dtmax, param.dt);
   if ( mpi.Root() )
   {
      std::cout << "Number of Time Steps:  " << nsteps << std::endl;
      std::cout << "Time Step Size:        " << param.dt / tScale_ << "ns" << std::endl;
   }

   // Create the ODE solver
   EHSVSolver siaSolver(param.integrator_order, VlasovMaxwell,VlasovMaxwell.GetNegCurl(), VlasovMaxwell);

   // Initialize VisIt visualization
   mfem::VisItDataCollection visit_dc("MaxP-Parallel", &pmesh);
   // Initialize ParaView visualization
   mfem::ParaViewDataCollection paraview_dc("MaxP-Parallel", &pmesh);

   double t = param.time_ini;
   VlasovMaxwell.SetTime(t);

   if ( param.visit )
   {
      Fields.RegisterVisItFields(visit_dc);
      Fields.WriteVisItFields(0);
   }

   if ( param.paraview )
   {
      Fields.RegisterParaViewFields(paraview_dc);
      Fields.WriteParaViewFields(0);

   }

   // The main time evolution loop.
   int it = 1;
   while (t < param.time_fin)
   {
      // Run the simulation until a snapshot is needed
      siaSolver.Run(Fields.GetBField(), Fields.GetEField(),
         t, param.dt, std::max(t + param.dt, param.time_ini + param.time_diag * it));

      // Approximate the current energy if the fields
      energy = VlasovMaxwell.GetEnergy();
      if ( mpi.Root() )
      {
         std::cout << "Energy(" << t/tScale_ << "ns):  " << energy << "J" << std::endl;
      }

      // Update local DoFs with current true DoFs
      Fields.SyncGridFuncs();

      // Write fields to disk for VisIt
      if ( param.visit )
      {
         Fields.WriteVisItFields(it, t);
      }

      // Write fields to disk for ParaView
      if ( param.paraview )
      {
         Fields.WriteParaViewFields(it, t);
      }

      it++;
   }

   return 0;
}

int
SnapTimeStep(double tmax, double dtmax, double & dt)
{
   double dsteps = tmax/dtmax;

   int nsteps = pow(10,(int)ceil(log10(dsteps)));

   for (int i=1; i<=5; i++)
   {
      int a = (int)ceil(log10(dsteps/pow(5.0,i)));
      int nstepsi = (int)pow(5,i)*std::max(1,(int)pow(10,a));

      nsteps = std::min(nsteps,nstepsi);
   }

   dt = tmax / nsteps;

   return nsteps;
}
