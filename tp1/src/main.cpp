/*---------------------------------------------------------------------------*/
/* "main.cpp"                                                                */
/*                                                                           */
/* HyPERION entry point.                                                     */
/*---------------------------------------------------------------------------*/

#include <iostream>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <yaml-cpp/yaml.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "HyperionMainDriver.hpp"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
  if (argc < 1) {
    std::cout << "Usage : h2p dataset.yaml\n";
    return 1;
  }

  // Read the dataset file
  std::string dataset_file(argv[1]);

  YAML::Node dataset;
  try {
    dataset = YAML::LoadFile(dataset_file);
  } catch (const YAML::Exception &e) {
    throw std::runtime_error("Error: can't load file " + dataset_file);
  }

  HyperionMainDriver driver(dataset);
  driver.load_mesh();
  driver.run();

  return 0;
}
