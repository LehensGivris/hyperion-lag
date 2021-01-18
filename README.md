# HyPERION-lag

This is **HyPERION**, as in *Hydrodynamics Platform for Exascale Research, In situ analysis and OptimizatioN* - Lagrangian variant

The example used in this project is heavily inspired by
[MicroHydroProject](https://github.com/cea-hpc/Modane/tree/master/plugins/fr.cea.modane.ui/examples/MicroHydroProject)
in the [Modane](https://github.com/cea-hpc/Modane) project from CEA.

## Compilation

A compiler supporting the C++ 17 standard is required. A recent CMake version is recommended. As for third-party
libraries, GMSH and VTK are both required. Additionally, ParaView compiled with Catalyst is needed for In Situ runs.

    cmake -S . -B build -DGMSH_DIR=/path/to/gmsh -DVTK_DIR=/path/to/vtk
    cmake --build build

To enable In Situ analysis, add the following options:

    cmake -S . -B build [...] -DParaView_DIR=/path/to/paraview -DHYPERION_ENABLE_INSITU=ON

## Usage

In the test case directory, for example `test/sod`, run :

    /path/to/h2p sod2d.yaml

VTK output files are produced in the current directory.

For In Situ runs, a python script describing the analysis pipeline is required.
As an example, `catalyst_insitu.py` has been created in `test/sod`.
