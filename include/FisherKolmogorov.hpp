#ifndef FISCHER_KOLMOGOROV_HPP
#define FISCHER_KOLMOGOROV_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>
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

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <boost/geometry.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "BrainRegions.hpp"

using namespace dealii;

/**
 * Class representing the Fisher-Kolmogorov equation.
 */
class FisherKolmogorov
{
public:
  static constexpr unsigned int dim = 3;

  /**
   * Function for the isotropic diffusion coefficient dext for
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
   * Function for the isotropic diffusion coefficient dext for
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
   * Function for the anisotropic diffusion coefficient daxn for
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
   * Function for the anisotropic diffusion coefficient daxn for
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
   * Function for the growth coefficient alpha for white matter.
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
   * Function for the growth coefficient alpha for grey matter.
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
   * Function for the initial concentration of the misfolded
   * protein c(t = 0), depending on the seeding region of the problem.
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

  // Number of MPI processes.
  const unsigned int mpi_size;

  // Rank of the current MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Coefficient dext for white matter.
  IsotropicDiffusionCoefficientWhite isotropic_coefficient_white;

  // Coefficient dext for grey matter.
  IsotropicDiffusionCoefficientGrey isotropic_coefficient_grey;

  // Coefficient daxn for white matter.
  AnisotropicDiffusionCoefficientGrey anisotropic_coefficient_grey;

  // Coefficient daxn for grey matter.
  AnisotropicDiffusionCoefficientWhite anisotropic_coefficient_white;

  // Coefficient alpha for white matter.
  GrowthCoefficientWhite growth_coefficient_white;

  // Coefficient alpha for grey matter.
  GrowthCoefficientGrey growth_coefficient_grey;

  // Initial concentration c(t = 0).
  FunctionC0 c_0;

  // Fiber orientation
  std::string orientation;

  // Time t.
  double time;

  // Final time.
  const double T;

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Time step.
  const double deltat;

  // Center of the domain.
  Point<dim> global_center;

  // Grey matter region.
  Grey_matter<dim> grey_matter;

  // Direction function of the fiber orientation.
  std::unique_ptr<Function<dim>> direction;

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
