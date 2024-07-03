### Organizing the source code
Please place all your sources into the `src` folder.

Binary files must not be uploaded to the repository (including executables).

Mesh files should not be uploaded to the repository. If applicable, upload `gmsh` scripts with suitable instructions to generate the meshes (and ideally a Makefile that runs those instructions). If not applicable, consider uploading the meshes to a different file sharing service, and providing a download link as part of the building and running instructions.

## Build with CMake
This project uses CMake as its build system, follow the steps below to build the project using it.

### Prerequisites

- CMake (version 3.12 or higher)
- MPI
- deal.II (version 9.3.1)
- Boost (1.72.0)

### Building the project

To build the executable, make sure you have loaded the needed modules with
```bash
$ module load gcc-glibc dealii
```
Then run the following commands:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```
Two executables are created into `build`, 'main1D' for 1D simulation and 'main' for 2D and 3D simulations. 

### Mesh Generation
Gmsh is required in order to generate the meshes (mesh1D.msh, ellipse.msh, brain-h3.0.msh) starting from the .geo files located in the mesh folder.

Use the Makefile running the commands:
```bash
$ cd mesh
$ make
```

### Run the executables

'main1D' can be executed through
```bash
$ ./main1D
```
Since 'main' supports parallel execution with MPI paradigm via deal.ii, it is possible to execute it in two ways:
1. **Serial Version**
```bash
$ ./main
```
2. **Parallel Version** (make sure you have enough hardware resources to run with the number of threads you entered)
```bash
$ mpirun -n <#threads> ./main
```

## Output
It is suggested to view the output files via Paraview, which are created after execution within the `build` folder.

'main1D' creates file with .vtk extension, instead 'main' creates .vtu files and a .pvtu file (open the .pvtu file on Paraview).



