# Finite element solver for the Fisher Kolmogorov equation
This project is a finite elements solver designed to numerically solve the Fisher-Kolmogorov equation in the domain of neurodegenerative-diseases.
The equation solved is:

$$
\begin{cases}
    \frac{\partial c}{\partial t} - \nabla \cdot (\mathbf{D} \nabla c) - \alpha c(1 - c) = 0 & \text{in } \Omega \times (0, T], \\
    \mathbf{D} \nabla c \cdot \mathbf{n} = 0 & \text{on } \partial \Omega \times (0, T], \\ 
    c(t = 0) = c_0 & \text{in } \Omega \times \{0\}.
\end{cases}
$$

with domain $\Omega \subset \mathbb{R}^d$, spatial coordinates $\textbf{x} \in \Omega^d$ and where 
- $c(\textbf{x},t): \Omega \times (0,+\infty)  \rightarrow \mathbb{R}$ represents the relative concentration of misfolded protein; 
- $\alpha \in \mathbb{R}$  characterizes the growth of the concentration of misfolded protein; 
- $\textbf{D} \ \in \mathbb{R}^{d\times d}$ represents the spreading of misfolded protein; 
- ${\textbf{n}}$ $\in \mathbb{R}^d$ is the axonal fiber direction.


To model the two different spreading mechanisms of misfolded proteins across the brain, the tensor $\textbf{D}$ is defined as: $\textbf{D} = d^{ext} \textbf{I} + d^{axn} \textbf{n} \otimes \textbf{n} $
where the extracellular diffusion $d^{ext}$ is associated with the isotropic diffusion of misfolded protein through the extracellular space, while the axonal transport $d^{axn}$ is associate with anisotropic diffusion of misfolded protein along the local axonal direction $\textbf{n}$.

There are two versions of the solver: a 1D version, used for convergence studies, and a general solver that accounts for different dimensions of the problem and offers more realistic physical modeling. We recommend working with the latter.


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

## Usage
### Set parameters  


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



