/*---------------------------------------------------------------------------*/
/* "HyperionMainDriver.hpp"                                                  */
/*                                                                           */
/* HyPERION driver.                                                          */
/*---------------------------------------------------------------------------*/
#ifndef HYPERION_HYPERIONMAINDRIVER_HPP
#define HYPERION_HYPERIONMAINDRIVER_HPP
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <vtkSmartPointer.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class vtkUnstructuredGrid;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class HyperionMainDriver
{
  // Boundary conditions groups dimension
  inline static constexpr int BC_ENV_DIM = 1;
  // Initial conditions groups dimension
  inline static constexpr int IC_ENV_DIM = 2;
  // GMSH quad cell id
  inline static constexpr int MSH_QUAD_4 = 3;

public:
  HyperionMainDriver(const YAML::Node dataset);
  ~HyperionMainDriver() = default;

public:
  void load_mesh();
  int run();

private:
  /// Dataset root node
  YAML::Node m_dataset;

  /// Computational grid
  vtkSmartPointer<vtkUnstructuredGrid> m_mesh;

  /// Environnment for cell
  std::map<int, std::string> m_cell_envs;
  /// Environments for node
  std::map<int, std::vector<std::string>> m_node_envs;

  /// Lookup tables for cell ids
  std::map<int, int> m_msh_vtk_cells;
  std::map<int, int> m_vtk_msh_cells;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif //HYPERION_HYPERIONMAINDRIVER_HPP
