/*---------------------------------------------------------------------------*/
/* "HydroVars.hpp"                                                           */
/*                                                                           */
/* Hydro variables.                                                          */
/*---------------------------------------------------------------------------*/

#include <iostream>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <yaml-cpp/yaml.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "HydroVars.hpp"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

HydroVars::HydroVars(int nb_cells, int nb_nodes)
  : m_nb_cells(nb_cells), m_nb_nodes(nb_nodes)
{
  this->init(nb_cells, nb_nodes, 0.0, 0.0,
             0.0, 0.0, 1.4, 0.0, 0.0, 0.0, 0.0,
             std::make_pair(0.0, 0.0), std::make_pair(0.0, 0.0));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void HydroVars::init(int nb_cells, int nb_nodes,
                     double density, double pressure,
                     double internal_energy,
                     double sound_speed,double adiabatic_cst,
                     double artificial_viscosity,
                     double cell_volume, double cell_mass,
                     double node_mass,
                     std::pair<double, double> velocity,
                     std::pair<double, double> force)
{
  std::fill_n(std::back_inserter(m_density), nb_cells, density);
  std::fill_n(std::back_inserter(m_old_density), nb_cells, density);
  std::fill_n(std::back_inserter(m_pressure), nb_cells, pressure);
  std::fill_n(std::back_inserter(m_internal_energy), nb_cells, internal_energy);
  std::fill_n(std::back_inserter(m_sound_speed), nb_cells, sound_speed);
  std::fill_n(std::back_inserter(m_cell_volume), nb_cells, cell_volume);
  std::fill_n(std::back_inserter(m_old_cell_volume), nb_cells, cell_volume);
  std::fill_n(std::back_inserter(m_cell_mass), nb_cells, cell_mass);
  std::fill_n(std::back_inserter(m_adiabatic_cst), nb_cells, adiabatic_cst);
  std::fill_n(std::back_inserter(m_artificial_viscosity), nb_cells, artificial_viscosity);

  std::fill_n(std::back_inserter(m_node_mass), nb_nodes, node_mass);
  std::fill_n(std::back_inserter(m_velocity), nb_nodes, velocity);
  std::fill_n(std::back_inserter(m_force), nb_nodes, force);

  std::cout << "[HydroVars::init] Initialized variables\n";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void HydroVars::setup_sod(const YAML::Node dataset,
                          std::map<int, std::string> cell_envs,
                          std::map<int, int> vtk_msh_cells)
{
  std::map<std::string, YAML::Node> initial_conditions;
  YAML::Node ics = dataset["InitialConditions"];
  for (YAML::const_iterator it = ics.begin(); it != ics.end(); ++it) {
    const YAML::Node& ic = *it;
    initial_conditions[ic["Name"].as<std::string>()] = ic;
  }

  for (int c = 0; c < m_nb_cells; ++c) {
    std::string env_name = cell_envs[vtk_msh_cells[c]];
    YAML::Node env_ic = initial_conditions[env_name];
    m_pressure[c] = env_ic["Pressure"].as<double>();
    m_density[c] = env_ic["Density"].as<double>();
    m_old_density[c] = env_ic["Density"].as<double>();
  }

  std::cout << "[HydroVars::setup_sod] Sod problem setup\n";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
