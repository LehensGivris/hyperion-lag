/*---------------------------------------------------------------------------*/
/* SPDX-License-Identifier: GPL-3.0-only                                     */
/* Copyright (C) 2020 Benjamin Fovet                                         */
/*---------------------------------------------------------------------------*/
/* "CatalystInsituAdaptor.hpp"                                               */
/*                                                                           */
/* TODO : description.                                                       */
/*---------------------------------------------------------------------------*/

#include <iostream>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkUnstructuredGrid.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "CatalystInsituAdaptor.hpp"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace CatalystAdaptor
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

vtkCPProcessor *processor = nullptr;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void initialize(const std::string& script)
{
  if (processor == nullptr) {
    processor = vtkCPProcessor::New();
    processor->Initialize();
  } else {
    processor->RemoveAllPipelines();
  }

  vtkNew<vtkCPPythonScriptPipeline> pipeline;
  pipeline->Initialize(script.c_str());
  processor->AddPipeline(pipeline);

  std::cout << "[CatalystInsituAdaptor::init] Initialized Catalyst In Situ Adaptor\n";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void finalize()
{
  if (processor) {
    processor->Delete();
    processor = nullptr;
  }

  std::cout << "[CatalystInsituAdaptor::finalize] Finalized Catalyst In Situ Adaptor\n";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void execute(vtkSmartPointer<vtkUnstructuredGrid> grid, double time,
             int timestep, bool last_timestep)
{
  vtkNew<vtkCPDataDescription> data_description;
  data_description->AddInput("input");
  data_description->SetTimeData(time, timestep);

  if (last_timestep) {
    data_description->ForceOutputOn();
  }

  if (processor->RequestDataDescription(data_description) != 0) {
    vtkCPInputDataDescription* idd = data_description->GetInputDescriptionByName("input");
    idd->SetGrid(grid);
    processor->CoProcess(data_description);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
