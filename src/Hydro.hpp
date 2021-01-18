/*---------------------------------------------------------------------------*/
/* "Hydro.hpp"                                                               */
/*                                                                           */
/* Hydrodynamics methods.                                                    */
/*---------------------------------------------------------------------------*/
#ifndef HYPERION_HYDRO_HPP
#define HYPERION_HYDRO_HPP
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <vtkSmartPointer.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class HydroVars;
class vtkUnstructuredGrid;
class vtkXMLUnstructuredGridWriter;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class Hydro
{
public:
  Hydro(YAML::Node dataset,
        vtkSmartPointer<vtkUnstructuredGrid> mesh,
        HydroVars* vars);
  ~Hydro() = default;

public:
  void init();
  void finalize();

  void compute_volume();
  void compute_cqs(std::pair<double, double>* coord,
                   int cell_idx, int nb_nodes_of_cell);
  void compute_pressure_force();
  void compute_artificial_viscosity();
  void compute_velocity();
  void apply_boundary_condition(const std::map<int, std::vector<std::string>>& node_envs);
  void move_nodes();
  void update_density();
  void apply_eos();
  void compute_dt();

  void dump(int step, double simulation_time);
  void analyze_insitu(double simulation_time, int iteration, bool last_iteration);

public:
  void set_dt(double dt) { m_dt = dt; }
  double dt() { return m_dt; }

private:
  void update_fields(double simulation_time);
  void init_insitu();
  void finalize_insitu();

private:
  YAML::Node m_dataset;
  vtkSmartPointer<vtkUnstructuredGrid> m_mesh;
  HydroVars* m_vars;
  vtkNew<vtkXMLUnstructuredGridWriter> m_writer;

  double m_dt;
  double m_dt_staggered;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif //HYPERION_HYDRO_HPP
