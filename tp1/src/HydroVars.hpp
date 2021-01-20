/*---------------------------------------------------------------------------*/
/* "HydroVars.hpp"                                                           */
/*                                                                           */
/* Hydro variables.                                                          */
/*---------------------------------------------------------------------------*/
#ifndef HYPERION_HYDROVARS_HPP
#define HYPERION_HYDROVARS_HPP
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <vector>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

struct HydroVars
{
  HydroVars(int nb_cells, int nb_nodes);
  void init(int nb_cells, int nb_nodes,
            double density, double pressure, double internal_energy,
            double sound_speed, double cell_volume, double cell_mass,
            double adiabatic_cst, double artificial_viscosity, double node_mass,
            std::pair<double, double> velocity, std::pair<double, double> force);

  void setup_sod(const YAML::Node dataset,
                 std::map<int, std::string> cell_envs,
                 std::map<int, int> vtk_msh_cells);

  // Cell centered scalar variables
  std::vector<double> m_density;
  std::vector<double> m_old_density;
  std::vector<double> m_pressure;
  std::vector<double> m_internal_energy;
  std::vector<double> m_sound_speed;
  std::vector<double> m_adiabatic_cst;
  std::vector<double> m_artificial_viscosity;

  std::vector<double> m_cell_volume;
  std::vector<double> m_old_cell_volume;
  std::vector<double> m_cell_mass;

  // Node centered scalar variables
  std::vector<double> m_node_mass;

  // Vector field variables
  std::vector<std::pair<double, double>> m_node_coord;
  std::vector<std::pair<double, double>> m_velocity;
  std::vector<std::pair<double, double>> m_force;
  std::vector<std::vector<std::pair<double, double>>> m_cqs;

  // Mesh properties
  int m_nb_cells;
  int m_nb_nodes;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif //HYPERION_HYDROVARS_HPP
