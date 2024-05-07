#ifndef FISCHER_KOLMOGOROV_HPP
#define FISCHER_KOLMOGOROV_HPP

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

// Two header taken by reading documentation in order to obatin finite element
// object
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

// to deal with linear systems (GMRES solver in combination with AMG
// preconditioner)
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

// Header taken in order to obatin ConditionalOStream type
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/numerics/data_out.h>     //used in output method
#include <deal.II/numerics/vector_tools.h> //used in solve method

#include <filesystem>
#include <fstream>
#include <iostream>


using namespace dealii;

#include <string>

/*
// SeedingRegion is a base class for different seeding regions
class SeedingRegion
{
public:

  SeedingRegion(double x_min_,
                double x_max_,
                double y_min_,
                double y_max_,
                double z_min_,
                double z_max_)
    : x_min(x_min_)
    , x_max(x_max_)
    , y_min(y_min_)
    , y_max(y_max_)
    , z_min(z_min_)
    , z_max(z_max_)
  {}

  protected:
  double x_min, x_max;
  double y_min, y_max;
  double z_min, z_max;
};

// todo
class Tau_inclusions : public SeedingRegion
{
public:
  Tau_inclusions()
    : SeedingRegion(63.0, 80.0, 48.0, 60.0, 50.0, 67.0)
  {}
};
// todo
class Amyloid_Beta_deposits : public SeedingRegion
{
public:
  Amyloid_Beta_deposits()
    : SeedingRegion(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  {}
};

// todo
class TPD43_inclusions : public SeedingRegion
{
public:
  TPD43_inclusions()
    : SeedingRegion(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  {}
};

// Factory function to create seeding regions
SeedingRegion *
getseddingRegion(std::string region)
{
  if (region == "Tau inclusions")
    return new Tau_inclusions();
  else if (region == "Amyloid-Beta deposits")
    return new Amyloid_Beta_deposits();
  else if (region == "TPD-43 inclusions")
    return new TPD43_inclusions();
  else
    throw std::invalid_argument("Invalid seeding region");
}
*/

class FisherKolmogorov
{
public:
  // Physical dimension (1D, 3D)
  static constexpr unsigned int dim = 3;
  class GrowthCoefficient : public Function<dim>
  {
  public:
    GrowthCoefficient()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.0;
    }
  };

  class FunctionC0 : public Function<dim>
  {
  public:
    FunctionC0()
    {}

    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      // tau inclusions
      if ((p[0] > 63 && p[0] < 80) && (p[1] > 62 && p[1] < 90) &&
          (p[2] > 46 && p[2] < 67))
        return 0.1;
      else
        return 0.0;
    }
  };

  FisherKolmogorov(const std::string  &mesh_file_name_,
                   const unsigned int &r_,
                   const double       &deltat_,
                   const double       &T_,
                   const double       &dext_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , T(T_)
    , mesh_file_name(mesh_file_name_)
    , r(r_)
    , deltat(deltat_)
    , mesh(MPI_COMM_WORLD)
  {
    spreading_coefficient.clear();
    for (size_t i = 0; i < dim; i++)
      spreading_coefficient[i][i] = dext_;
    // method to choose the region at runtime based on the name passed.
  }

  // setup the problem
  void
  setup();

  // Solve the problem.
  void
  solve();


protected:
  // Assemble the tangent problem.
  void
  assemble_system();

  // Solve the linear system associated to the tangent problem.
  void
  solve_linear_system();

  // Solve the problem for one time step using Newton's method.
  void
  solve_newton();

  // Output
  void
  output(const unsigned int &time_step) const;

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Coefficient alpha
  GrowthCoefficient growth_coefficient;

  // spreading coefficient D
  Tensor<2, dim> spreading_coefficient;

  // Initial concentration c(t = 0)
  FunctionC0 c_0;

  double time;

  const double T;


  // SeedingRegion *region = nullptr;

  // Path to the mesh file.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  const double deltat;

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

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // Jacobian matrix.
  TrilinosWrappers::SparseMatrix jacobian_matrix;

  // Residual vector.
  TrilinosWrappers::MPI::Vector residual_vector;

  // Increment of the solution between Newton iterations.
  TrilinosWrappers::MPI::Vector delta_owned;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::Vector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::Vector solution;

  // System solution at previous time step.
  TrilinosWrappers::MPI::Vector solution_old;
};

#endif
