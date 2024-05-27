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

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

// Header taken in order to obatin ConditionalOStream type
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/numerics/data_out.h>     //used in output method
#include <deal.II/numerics/vector_tools.h> //used in solve method

#include <boost/geometry.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "BrainRegions.hpp"

using namespace dealii;

/**
 * This class is used to define the Fisher-Kolmogorov equation and solve it
 * numerically.
 */
class FisherKolmogorov
{
public:
  static constexpr unsigned int dim = 3;

  /**
   * This class is used to define the spreading coefficient D.
   */
  class SpreadingCoefficient
  {
  public:
    SpreadingCoefficient()
    {}

    /**
     * Compute the spreading coefficient D as:
     * \[
     * D = d_{ext}} \cdot \mathbf{I} + d_{axn} \cdot \mathbf{n} \otimes
     * \mathbf{n}
     * \]
     *
     * @param dext isotropic diffusion coefficient
     * @param daxn anisotropic diffusion coefficient
     * @param direction fiber orientation
     * @return Tensor<2, dim> the spreading coefficient D
     */
    Tensor<2, dim>
    value(const double         &dext,
          const double         &daxn,
          const Tensor<1, dim> &direction) const
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
  };

  /**
   * This class is used to define the isotropic diffusion coefficient dext for
   * white matter.
   */
  class IsotropicDiffusionCoefficientWhite : public Function<dim>
  {
  public:
    IsotropicDiffusionCoefficientWhite()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.5;
    }
  };

  /**
   * This class is used to define the isotropic diffusion coefficient dext for
   * grey matter.
   */
  class IsotropicDiffusionCoefficientGrey : public Function<dim>
  {
  public:
    IsotropicDiffusionCoefficientGrey()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.5;
    }
  };

  /**
   * This class is used to define the anisotropic diffusion coefficient daxn for
   * white matter.
   */
  class AnisotropicDiffusionCoefficientWhite : public Function<dim>
  {
  public:
    AnisotropicDiffusionCoefficientWhite()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 3.0;
    }
  };

  /**
   * This class is used to define the anisotropic diffusion coefficient daxn for
   * grey matter.
   */
  class AnisotropicDiffusionCoefficientGrey : public Function<dim>
  {
  public:
    AnisotropicDiffusionCoefficientGrey()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  /**
   * This class is used to define the growth coefficient alpha for white matter.
   */
  class GrowthCoefficientWhite : public Function<dim>
  {
  public:
    GrowthCoefficientWhite()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.6;
    }
  };

  /**
   * This class is used to define the growth coefficient alpha for grey matter.
   */
  class GrowthCoefficientGrey : public Function<dim>
  {
  public:
    GrowthCoefficientGrey()
    {}

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.3;
    }
  };

  /**
   * This class is used to define the initial concentration c(t = 0).
   */
  class FunctionC0 : public Function<dim>
  {
  public:
    FunctionC0(const std::string &seeding_region_)
      : seeding_region(std::move(getSeedingRegion<dim>(seeding_region_)))
    {}

    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      if (seeding_region->check_region(p))
        {
          return 0.1;
        }
      else
        return 0.0;
    }

  private:
    const std::unique_ptr<SeedingRegion<dim>> seeding_region;
  };

  FisherKolmogorov(const std::string  &mesh_file_name_,
                   const unsigned int &r_,
                   const double       &deltat_,
                   const double       &T_,
                   const std::string  &seeding_region_,
                   const std::string  &orientation_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , c_0(seeding_region_)
    , orientation(orientation_)
    , T(T_)
    , mesh_file_name(mesh_file_name_)
    , r(r_)
    , deltat(deltat_)
    , mesh(MPI_COMM_WORLD)
  {}

  void
  setup();

  void
  solve();

protected:
  void
  assemble_system();

  void
  solve_linear_system();

  void
  solve_newton();

  void
  output(const unsigned int &time_step) const;

  const unsigned int mpi_size;

  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  GrowthCoefficientWhite growth_coefficient_white;

  GrowthCoefficientGrey growth_coefficient_grey;

  SpreadingCoefficient spreading_coefficient;

  IsotropicDiffusionCoefficientWhite isotropic_coefficient_white;

  IsotropicDiffusionCoefficientGrey isotropic_coefficient_grey;

  AnisotropicDiffusionCoefficientGrey anisotropic_coefficient_grey;

  AnisotropicDiffusionCoefficientWhite anisotropic_coefficient_white;

  FunctionC0 c_0;

  // Fiber orientation
  std::string orientation;

  double time;

  const double T;

  const std::string mesh_file_name;

  const unsigned int r;

  const double deltat;

  Point<dim> global_center;

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
