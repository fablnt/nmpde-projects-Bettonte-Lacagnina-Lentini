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

There are two versions of the solver: a 1D serial version, used for convergence studies, and a general solver that accounts for different dimensions of the problem and offers more realistic physical modeling and computing performances. We recommend working with the latter.


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
Two executables are created into `build`: `main1D` for 1D simulation and `main` for 2D and 3D simulations. 

### Mesh Generation
Gmsh is required in order to generate the meshes (mesh1D.msh, ellipse.msh, brain-h3.0.msh) starting from the .geo files located in the `mesh` folder.

Use the Makefile running the commands:
```bash
$ cd mesh
$ make
```

## Usage
### Set parameters  
The parameters for the solver are set directly in the code. To modify them you need to modify the following fields:

- The following parameters can be modified directly in the `main.cpp` file:

| Parameter        | Description                                           | Example Value                |
|------------------|-------------------------------------------------------|------------------------------|
| `mesh_filename`  | Path to the mesh file                                 | `"../mesh/brain-h3.0.msh"`   |
| `degree`         | Degree of the finite element space                    | `1`                          |
| `deltat`         | Time step size                                        | `0.4`                        |
| `T`              | Final time for simulation                             | `40`                         |
| `seeding_region` | Region where protein seeding begins                   | `"Amyloid-Beta deposits"`    |
| `orientation`    | Fiber orientation                                     | `"axon-based"`               |

The available values for the seeding region variable are : "Tau inclusions", "TPD-43 inclusions", "Amyloid-Beta deposits". The available values for the fiber orientation are: "radial", "circumferential", "axon-based".

- To set the dimension of the problem you need to modify the following variable ([Line 54](include/FisherKolmogorov.hpp#L54)) in the `FisherKolmogorov.hpp` file:
```cpp
static constexpr unsigned int dim = 3;
```

- To set the coefficients of the Fisher-Kolmogorov equation, you need to modify the return value of the following method:
```cpp
 virtual double
    value(const Point<dim> & p,
          const unsigned int component = 0) const override
```
that is overridden in the classes representing the equation coefficients defined in the `FisherKolmogorov.hpp` file. In particular:
| Class                                                                            | Description                                                  |
|----------------------------------------------------------------------------------|--------------------------------------------------------------|
| `IsotropicDiffusionCoefficientWhite`                                             | Coefficient $d^{ext}$ for white matter region [(Lines 60-72)](include/FisherKolmogorov.hpp#L60-L72)      |
| `IsotropicDiffusionCoefficientGrey`                                              | Coefficient $d^{ext}$ for grey matter region [(Lines 78-90)](include/FisherKolmogorov.hpp#L78-L90)       |
| `AnisotropicDiffusionCoefficientWhite`                                           | Coefficient $d^{axn}$ for white matter region [(Lines 96-108)](include/FisherKolmogorov.hpp#L96-L108)    |
| `AnisotropicDiffusionCoefficientGrey`                                            | Coefficient $d^{axn}$ for grey matter region [(Lines 114-126)](include/FisherKolmogorov.hpp#L114-L126)   |
| `GrowthCoefficientWhite`                                                         | Coefficient $\alpha$ for white matter region [(Lines 131-143)](include/FisherKolmogorov.hpp#L131-L143)   |
| `GrowthCoefficientGrey`                                                          | Coefficient $\alpha$ for grey matter region [(Lines 148-160)](include/FisherKolmogorov.hpp#L148-L160)    |
| `FunctionC0`                                                                     | Initial concentration $c(t = 0)$ [(Lines 166-187)](include/FisherKolmogorov.hpp#L166-L187)               |



### Run the executables

`main1D` can be executed through
```bash
$ ./main1D
```
Since `main` supports parallel execution with MPI paradigm via deal.ii, it is possible to execute it in two ways:
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



