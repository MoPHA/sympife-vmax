
#include "parameters.hpp"


using namespace std;

ParametersInput::ParametersInput(const char *config_file)
: N_curlB(1.0),
N_curlE(1.0),
N_eta(0.0),
time_ini(0.0),
time_fin(1.0),
time_diag(1.0),
dt(0.1),
dt_safe(0.95),
epsilon0(8.8541878176e-12),
mu0(4.0e-7*M_PI),
rho0(1.),
fe_order(1),
integrator_order(1),
mesh_serial_ref(1),
mesh_para_ref(1),
num_species(1),
case_name(""),
mesh_path(""),
visit(false),
paraview(false)
{
  // std::ifstream is RAII, i.e. no need to call close

  ifstream cFile(config_file);
  if (cFile.is_open())
  {
    string line;
    while(getline(cFile, line)){
      line.erase(remove_if(line.begin(), line.end(), ::isspace),
      line.end());
      if(line[0] == '#' || line.empty())
      continue;
      auto delimiterPos = line.find("=");
      string name = line.substr(0, delimiterPos);
      auto value = line.substr(delimiterPos + 1);
      replace(value.begin(),value.end(), ',', ' ');
      stringstream value_stream(value);
      int itemp;
      uint uitemp;
      double dtemp;
      bool btemp;

      //double parameters
      if (!name.find("N_curlB"))
      {
        value_stream >> this->N_curlB;
      }
      else if (!name.find("N_curlE"))
      {
        value_stream >> this->N_curlE;
      }
      else if (!name.find("N_eta"))
      {
        value_stream >> this->N_eta;
      }
      else if (!name.find("time_ini"))
      {
        value_stream >> this->time_ini;
      }
      else if (!name.find("time_fin"))
      {
        value_stream >> this->time_fin;
      }
      else if (!name.find("time_diag"))
      {
        value_stream >> this->time_diag;
      }
      else if (!name.find("dt"))
      {
        value_stream >> this->dt;
      }
      else if (!name.find("dt_safe"))
      {
        value_stream >> this->dt_safe;
      }
      else if (!name.find("epsilon0"))
      {
        value_stream >> this->epsilon0;
      }
      else if (!name.find("mu0"))
      {
        value_stream >> this->mu0;
      }
      else if (!name.find("rho0"))
      {
        value_stream >> this->rho0;
      }
      //Array<double> parameters
      else if (!name.find("N_q_s"))
      {
        while (value_stream >> dtemp) this->N_q_s.Append(dtemp);
      }
      else if (!name.find("N_E_s"))
      {
        while (value_stream >> dtemp) this->N_E_s.Append(dtemp);
      }
      else if (!name.find("N_vxB_s"))
      {
        while (value_stream >> dtemp) this->N_vxB_s.Append(dtemp);
      }
      else if (!name.find("N_v_s"))
      {
        while (value_stream >> dtemp) this->N_v_s.Append(dtemp);
      }
      else if (!name.find("eps_sphere"))
      {
        while (value_stream >> dtemp) this->eps_sphere.Append(dtemp);
      }
      else if (!name.find("mu_shell"))
      {
        while (value_stream >> dtemp) this->mu_shell.Append(dtemp);
      }
      else if (!name.find("sigma_sphere"))
      {
        while (value_stream >> dtemp) this->sigma_sphere.Append(dtemp);
      }
      else if (!name.find("dipole_pulse"))
      {
        while (value_stream >> dtemp) this->dipole_pulse.Append(dtemp);
      }
      else if (!name.find("point_charges"))
      {
        while (value_stream >> dtemp) this->point_charges.Append(dtemp);
      }
      //uint parameters
      else if (!name.find("fe_order"))
      {
        value_stream >> this->fe_order;
      }
      else if (!name.find("integrator_order"))
      {
        value_stream >> this->integrator_order;
      }
      else if (!name.find("mesh_serial_ref"))
      {
        value_stream >> this->mesh_serial_ref;
      }
      else if (!name.find("mesh_para_ref"))
      {
        value_stream >> this->mesh_para_ref;
      }
      else if (!name.find("num_species"))
      {
        value_stream >> this->num_species;
      }
      //Array<int> parameters
      else if (!name.find("bc_absorbing"))
      {
        while (value_stream >> itemp) this->bc_absorbing.Append(itemp);
      }
      else if (!name.find("bc_dirichlet"))
      {
        while (value_stream >> itemp) this->bc_dirichlet.Append(itemp);
      }
      else if (!name.find("markers_numbers"))
      {
        while (value_stream >> itemp) this->markers_numbers.Append(itemp);
      }
      //string parameters
      else if (!name.find("case_name"))
      {
        value_stream >> this->case_name;
      }
      else if (!name.find("mesh_path"))
      {
        value_stream >> this->mesh_path;
      }
      //bool parameters
      else if (!name.find("visit"))
      {
        value_stream >> this->visit;
      }
      else if (!name.find("paraview"))
      {
        value_stream >> this->paraview;
      }
      else
      {
        cerr << "Parameter " << name << " is not recognised" << endl;
      }
    }
    num_species = N_v_s.Size();
  }
  else {
    cerr << "Couldn't open config file for reading."  << endl;
  }
}
