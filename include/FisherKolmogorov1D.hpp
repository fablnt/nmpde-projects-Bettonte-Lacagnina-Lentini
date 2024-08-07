#ifndef FISHER_KOLMOGOROV_1D_HPP
#define FISHER_KOLMOGOROV_1D_HPP

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
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * Class representing the Fisher-Kolmogorov equation in 1D, solved serially.
 */
class FisherKolmogorov1D
{
public:
  static constexpr unsigned int dim = 1;

  /**
   * Class representing the initial concentration c(t = 0).
   */
  class FunctionC0 : public Function<dim>
  {
  public:
    FunctionC0()
    {}

    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      if ((p[0] - 0) == 0)
        return 0.1;
      else
        return 0.0;
    }
  };

  FisherKolmogorov1D(const std::string  &mesh_file_name_,
                     const unsigned int &r_,
                     const double       &deltat_,
                     const double       &T_,
                     const double       &spreading_coeff_,
                     const double       &growth_coeff_)
    : growth_coefficient(growth_coeff_)
    , spreading_coefficient(spreading_coeff_)
    , c_0()
    , time(0.0)
    , T(T_)
    , mesh_file_name(mesh_file_name_)
    , r(r_)
    , deltat(deltat_)
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

  // Coefficient alpha
  const double growth_coefficient;

  // Coefficient D
  const double spreading_coefficient;

  // Initial concentration c(t = 0)
  FunctionC0 c_0;

  // Current time.
  double time;

  // Final time.
  const double T;

  // Path to the mesh file.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Time step.
  const double deltat;

  // Triangulation.
  Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // Sparsity pattern.
  SparsityPattern sparsity_pattern;

  // Jacobian matrix.
  SparseMatrix<double> jacobian_matrix;

  // Residual vector.
  Vector<double> residual_vector;

  // Increment of the solution between Newton iterations.
  Vector<double> delta;

  // System solution (including ghost elements).
  Vector<double> solution;

  // System solution at previous time step.
  Vector<double> solution_old;
};

#endif