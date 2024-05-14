#ifndef FISCHER_KOLMOGOROV_HPP
#define FISCHER_KOLMOGOROV_HPP

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_function.h>

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
#include <deal.II/base/tensor.h>

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
#include <memory>
#include <string>

#include "SeedingRegions.hpp"


using namespace dealii;


class FisherKolmogorov
{
public:
  // Physical dimension (1D, 3D)
  static constexpr unsigned int dim = 3;

  // Spreading coefficient D = dext * I + daxn * n x n, where x is a tensor
  // product and n is the fiber orientation
  class SpreadingCoefficient
  {
  public:
    SpreadingCoefficient(const double &dext_, const double &daxn_)
      : dext(dext_)
      , daxn(daxn_)
    {}

    // Returns the spreading coefficient D = dext * I + daxn * n x n
    Tensor<2, dim>
    value(const Point<dim> & /*p*/, const Tensor<1, dim> &direction) const
    {
      Tensor<2, dim> D;
      D.clear();
      for (unsigned int i = 0; i < dim; i++)
        D[i][i] = dext;

      Tensor<2, dim> S;
      S.clear();
      S = outer_product(direction, direction) * daxn;

      return D + S;
    }

  private:
    const double dext, daxn;
  };

  class GrowthCoefficientWhite : public Function<dim>
  {
  public:
    GrowthCoefficientWhite()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.5;
    }
  };

  class GrowthCoefficientGrey : public Function<dim>
  {
  public:
    GrowthCoefficientGrey()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.5;
    }
  };

  class FunctionC0 : public Function<dim>
  {
  public:
    FunctionC0(const std::string &seeding_region_)
      : seeding_region(std::move(getSeedingRegion(seeding_region_)))
    {}

    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      if ((p[0] > seeding_region->x_min && p[0] < seeding_region->x_max) &&
          (p[1] > seeding_region->y_min && p[1] < seeding_region->y_max) &&
          (p[2] > seeding_region->z_min && p[2] < seeding_region->z_max))
        {
          return 0.1;
        }
      else
        return 0.0;
    }

  private:
    const std::unique_ptr<SeedingRegion> seeding_region;
  };

  FisherKolmogorov(const std::string  &mesh_file_name_,
                   const unsigned int &r_,
                   const double       &deltat_,
                   const double       &T_,
                   const double       &dext_,
                   const double       &daxn_,
                   const std::string  &seeding_region_,
                   const std::string  &orientation_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , spreading_coefficient(dext_, daxn_)
    , c_0(seeding_region_)
    , orientation(orientation_)
    , T(T_)
    , mesh_file_name(mesh_file_name_)
    , r(r_)
    , deltat(deltat_)
    , mesh(MPI_COMM_WORLD)
  {}

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

  Tensor<1, dim>
  compute_radial_direction(
    const dealii::TriaActiveIterator<dealii::DoFCellAccessor<dim, dim, false>>
      &cell) const;

  Tensor<1, dim>
  compute_circumferential_direction(
    const dealii::TriaActiveIterator<dealii::DoFCellAccessor<dim, dim, false>>
      &cell) const;

  Tensor<1, dim>
  compute_axon_based_direction(
    const dealii::TriaActiveIterator<dealii::DoFCellAccessor<dim, dim, false>>
      &cell) const;

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Coefficient alpha
  GrowthCoefficientWhite growth_coefficient_white;
  GrowthCoefficientGrey  growth_coefficient_grey;

  // spreading coefficient D
  SpreadingCoefficient spreading_coefficient;

  // Initial concentration c(t = 0)
  FunctionC0 c_0;

  std::string orientation;

  double time;

  const double T;

  // Path to the mesh file.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  const double deltat;

  // Triangulation. The parallel::fullydistributed::Triangulation class manages
  // a triangulation that is completely distributed (i.e. each process only
  // knows about the elements it owns and its ghost elements).
  parallel::fullydistributed::Triangulation<dim> mesh;

  Point<dim> center;

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
