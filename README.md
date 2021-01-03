LehensGivris = Louis LE LANN (llelann) (CISD) Univ Bordeaux

Entraide avec Rémi LASVENES (rlasvenes) (CISD) Univ Bordeaux (Code pouvant être similaire)

# HyPERION-lag

This is **HyPERION**, as in *Hydrodynamics Platform for Exascale Research, In situ analysis and OptimizatioN* - Lagrangian variant

The example used in this project is heavily inspired by
[MicroHydroProject](https://github.com/cea-hpc/Modane/tree/master/plugins/fr.cea.modane.ui/examples/MicroHydroProject)
in the [Modane](https://github.com/cea-hpc/Modane) project from CEA.

## Compilation

A compiler supporting the C++ 17 standard is required. A recent CMake version is recommended. As for third-party
libraries, GMSH and VTK are both required.

    cmake -S . -B build -DGMSH_DIR=/path/to/gmsh -DVTK_DIR=/path/to/vtk
    cmake --build build

## Usage

In the test case directory, for example `test/sod`, run :

    /path/to/h2p sod2d.yaml

VTK output files are produced in the current directory.
