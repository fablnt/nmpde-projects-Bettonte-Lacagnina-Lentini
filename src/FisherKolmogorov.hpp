#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

//Two header taken by reading documentation in order to obatin finite element object
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

//Header taken in order to obatin ConditionalOStream type
#include <deal.II/base/conditional_ostream.h>

#include <fstream>
#include <iostream>
#include <filesystem>


using namespace dealii;

class FisherKolmogorov {

    public:

    // Physical dimension (1D, 3D)
    static constexpr unsigned int dim = 3;
    
    FisherKolmogorov(const std::string &mesh_file_name_, const unsigned int &r_)
    : mesh_file_name(mesh_file_name_)
    , r(r_)
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , mesh(MPI_COMM_WORLD)
    , pcout(std::cout, mpi_rank == 0)
  {}

   void setup();
   
   protected:
   // Path to the mesh file.
   const std::string mesh_file_name;
   
   // Polynomial degree.
   const unsigned int r;
   
   // Number of MPI processes.
   const unsigned int mpi_size;
   
   // This MPI process.
   const unsigned int mpi_rank;

  // Triangulation. The parallel::fullydistributed::Triangulation class manages
  // a triangulation that is completely distributed (i.e. each process only
  // knows about the elements it owns and its ghost elements).
  parallel::fullydistributed::Triangulation<dim> mesh;

    // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // Parallel output stream.
  ConditionalOStream pcout;
  
  // DoFs owned by current process.
  IndexSet locally_owned_dofs;


};