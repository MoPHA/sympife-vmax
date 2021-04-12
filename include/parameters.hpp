
#include "mfem.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>


#ifndef PARAM_IO
#define PARAM_IO


using namespace std;

class ParametersInput
{
public:

        ParametersInput(const char *config_file);

        double N_curlB,                 //Normalised curlB coefficient
        N_curlE,                        //Normalised curlE coefficient
        N_eta,                          //Normalised admittance for absorbing BC
        time_ini,                       //Initial time
        time_fin,                       //Final time
        time_diag,                      //Diagnostics time interval
        dt,                             //Time step (determined by GetMaximumTimeStep)
        dt_safe,                        //Safety scaling coefficent of the max time step
        epsilon0,                       //Vacuum permittivity
        mu0,                            //Vacuum permeability
        rho0;                             //Average density

        mfem::Array<double> N_q_s,      //Normalised current coefficients (1 per species)
        N_E_s,                          //Normalised electric force coefficients (1 per species)
        N_vxB_s,                        //Normalised magnetic force coefficients (1 per species)
        N_v_s,                          //Species normalised thermal speed
        eps_sphere,                     //Dielectric sphere parameters
        mu_shell,                       //Magnetic shell parameters
        sigma_sphere,                   //Conductive sphere parameters
        dipole_pulse,                   //Current dipole pulse parameters
        point_charges;                  //Point charges parameters

        uint fe_order,                  //Order of finite element basis functions (arbitrary)
        integrator_order,               //Order of symplectic time integrator (1 or even < 10)
        mesh_serial_ref,                //Number of refinement iterations on the serial mesh
        mesh_para_ref,                  //Number of refinement iterations on the decomposed mesh
        num_species;                    //Number of plasma species (inferred fron CJ, CVE and CVB)

        mfem::Array<int> bc_absorbing,  //Absorbing boundaries flag
        bc_dirichlet,                   //Dirichlet boundaries flag
        markers_numbers;                //Number of particles per species

        string case_name,               //Path to mesh file
        mesh_path;                      //Case ID for output labeling


        bool visit,
        paraview;
};


#endif
