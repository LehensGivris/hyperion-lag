/*---------------------------------------------------------------------------*/
/* "HyperionMainDriver.cpp"                                                  */
/*                                                                           */
/* HyPERION driver.                                                          */
/*---------------------------------------------------------------------------*/

#include <chrono>
#include <string>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <yaml-cpp/yaml.h>
#include <gmsh.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdList.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "HydroVars.hpp"
#include "Hydro.hpp"
#include "HyperionMainDriver.hpp"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

HyperionMainDriver::HyperionMainDriver(const YAML::Node dataset)
  : m_dataset(dataset)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void HyperionMainDriver::load_mesh()
{
  std::cout << "[Driver::load_mesh] Initialize GMSH API\n";
  gmsh::initialize();
  gmsh::option::setNumber("General.Terminal", 1);

  if (auto mesh_file_node = m_dataset["Mesh"]["MeshFile"]) {
    gmsh::open(mesh_file_node.as<std::string>());
  } else {
    throw std::runtime_error("[Driver::load_mesh] Error : 'MeshFile' property is required");
  }

  // Read environments for initial conditions
  std::vector<std::pair<int, int>> ic_envs;
  gmsh::model::getPhysicalGroups(ic_envs, IC_ENV_DIM);

  for (const auto& env : ic_envs) {
    int env_idx = env.second;
    std::string env_name;
    gmsh::model::getPhysicalName(2, env_idx, env_name);

    std::vector<int> entities;
    gmsh::model::getEntitiesForPhysicalGroup(IC_ENV_DIM, env_idx, entities);
    for (const auto& e : entities) {
      std::vector<std::size_t> cells;
      std::vector<std::size_t> nodes;
      gmsh::model::mesh::getElementsByType(MSH_QUAD_4, cells, nodes, e);
      for (const auto& c : cells) {
        m_cell_envs[c] = env_name;
      }
    }
  }

  // Read environments for boundary conditions
  std::vector<std::pair<int, int>> bc_envs;
  gmsh::model::getPhysicalGroups(bc_envs, BC_ENV_DIM);

  for (const auto& env : bc_envs) {
    int env_idx = env.second;
    std::string env_name;
    gmsh::model::getPhysicalName(BC_ENV_DIM, env_idx, env_name);

    std::vector<std::size_t> nodes;
    std::vector<double> coords;
    gmsh::model::mesh::getNodesForPhysicalGroup(BC_ENV_DIM, env_idx, nodes, coords);
    for (const auto& n : nodes) {
      m_node_envs[n].push_back(env_name);
    }
  }

  // Read node coordinates
  std::vector<std::size_t> nodes;
  std::vector<double> coords;
  [[maybe_unused]] std::vector<double> pcoords; // unused

  gmsh::model::mesh::getNodes(nodes, coords, pcoords);

  std::cout << "[Driver::load_mesh] Done reading mesh\n";
  std::cout << "[Driver::load_mesh] Initializing a VTK unstructured grid\n";

  // Create VTK points
  vtkNew<vtkPoints> points;
  points->SetDataTypeToDouble();

  // Insert points from Gmsh node coordinates
  for (std::size_t n = 0; n < nodes.size(); ++n) {
    points->InsertPoint(nodes[n] - 1, coords[n * 3], coords[n * 3 + 1], coords[n * 3 + 2]);
  }

  // Create a VTK unstructured grid
  m_mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_mesh->SetPoints(points);

  int nb_cells_to_allocate = 0;
  {
    std::vector<std::size_t> cells;
    std::vector<std::size_t> nodes;
    gmsh::model::mesh::getElementsByType(3, cells, nodes);
    nb_cells_to_allocate = cells.size();
  }

  // Allocate cells
  m_mesh->Allocate(nb_cells_to_allocate);

  // Get global cells and nodes
  nodes.clear();
  std::vector<std::size_t> cells;
  gmsh::model::mesh::getElementsByType(MSH_QUAD_4, cells, nodes);

  for (std::size_t c = 0; c < cells.size(); ++c) {
    m_msh_vtk_cells[cells[c]] = c;
    m_vtk_msh_cells[c] = cells[c];

    vtkNew<vtkIdList> cell_nodes;
    for (int n = 0; n < 4; ++n) {
      cell_nodes->InsertNextId(nodes[c * 4 + n] - 1);
    }
    m_mesh->InsertNextCell(VTK_QUAD, cell_nodes);
  }

  gmsh::finalize();

  std::cout << "[Driver::load_mesh] Mesh created\n";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int HyperionMainDriver::run()
{
  auto vars = new HydroVars(m_mesh->GetNumberOfCells(),
                            m_mesh->GetNumberOfPoints());
  vars->setup_sod(m_dataset, m_cell_envs, m_vtk_msh_cells);

  auto hydro = new Hydro(m_dataset, m_mesh, vars);
  hydro->init();

  std::cout << "[Driver::run] Simulation initialized, starting time loop\n";

  std::chrono::duration<double> computation_time(0.0);
  std::chrono::duration<double> iteration_time(0.0);
  double simulation_time = 0.0;
  int step = 0;

  double final_time = m_dataset["TimeManagement"]["FinalTime"].as<double>();
  std::cout << "[Driver::run] Final simulation time is " << final_time << " s\n";

  while (simulation_time <= final_time + hydro->dt()) {
    auto loop_start_time = std::chrono::high_resolution_clock::now();

    bool last_iteration = (simulation_time == final_time + hydro->dt());
    // TODO: if ANALYZE_INSITU, then call the appropriate method, else call dump

    hydro->compute_pressure_force();
    hydro->compute_artificial_viscosity();
    hydro->compute_velocity();
    hydro->apply_boundary_condition(m_node_envs);
    hydro->move_nodes();
    hydro->compute_volume();
    hydro->update_density();
    hydro->apply_eos();
    hydro->compute_dt();

    simulation_time += hydro->dt();
    step += 1;

    auto loop_end_time = std::chrono::high_resolution_clock::now();
    computation_time += loop_end_time - loop_start_time;
  }

  hydro->finalize();

  std::cout << "[Driver::run] Computation time : "
            << std::chrono::duration_cast<std::chrono::seconds>(computation_time).count()
            << " s\n";

  return 0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
