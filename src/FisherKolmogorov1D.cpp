#include "FisherKolmogorov1D.hpp"

void
FisherKolmogorov1D::setup()
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(mesh);

  std::ifstream grid_in_file(mesh_file_name);
  grid_in.read_msh(grid_in_file);

  std::cout << "  Number of elements = " << mesh.n_global_active_cells()
            << std::endl;

  // Initialize the finite element space.
  {
    std::cout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_Q<dim>>(r);

    std::cout << "  Degree                     = " << fe->degree << std::endl;
    std::cout << "  DoFs per cell              = " << fe->dofs_per_cell
              << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    std::cout << "  Quadrature points per cell = " << quadrature->size()
              << std::endl;
  }

  // Initialize the DoF handler.
  {
    std::cout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    std::cout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  // Initialize the linear system.
  {
    std::cout << "Initializing the linear system" << std::endl;

    std::cout << "  Initializing the sparsity pattern" << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    std::cout << "  Initializing the matrices" << std::endl;
    jacobian_matrix.reinit(sparsity_pattern);

    std::cout << "  Initializing the system right-hand side" << std::endl;
    residual_vector.reinit(dof_handler.n_dofs());
    std::cout << "  Initializing the solution vector" << std::endl;

    solution.reinit(dof_handler.n_dofs());
    solution_old = solution;

    delta.reinit(dof_handler.n_dofs());
  }
}

void
FisherKolmogorov1D::assemble_system()
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_residual(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  jacobian_matrix = 0.0;
  residual_vector = 0.0;

  // Value and gradient of the solution on current cell.
  std::vector<double>         solution_loc(n_q);
  std::vector<Tensor<1, dim>> solution_gradient_loc(n_q);

  // Value of the solution at previous timestep (un) on current cell.
  std::vector<double> solution_old_loc(n_q); // un+1

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell_matrix   = 0.0;
      cell_residual = 0.0;

      fe_values.get_function_values(solution, solution_loc);
      fe_values.get_function_gradients(solution, solution_gradient_loc);
      fe_values.get_function_values(solution_old, solution_old_loc);

      for (unsigned int q = 0; q < n_q; ++q)
        {
          // Evaluate coefficients on this quadrature node.
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Mass matrix.
                  cell_matrix(i, j) +=
                    fe_values.shape_value(i, q) * // time derivative term
                    fe_values.shape_value(j, q) / deltat * fe_values.JxW(q);

                  // Non-linear stiffness matrix, first term.
                  cell_matrix(i, j) +=
                    spreading_coefficient * fe_values.shape_grad(i, q) *
                    fe_values.shape_grad(j, q) * fe_values.JxW(q);

                  // Non-linear stiffness matrix, second term.
                  cell_matrix(i, j) -=
                    growth_coefficient * fe_values.shape_value(i, q) *
                    fe_values.shape_value(j, q) * fe_values.JxW(q);

                  // Non-linear stiffness matrix, third term.
                  cell_matrix(i, j) +=
                    2 * growth_coefficient * solution_loc[q] *
                    fe_values.shape_value(i, q) * fe_values.shape_value(j, q) *
                    fe_values.JxW(q);
                }

              // Assemble the residual vector (with changed sign).
              // Time derivative term.
              cell_residual(i) -= (solution_loc[q] - solution_old_loc[q]) /
                                  deltat * fe_values.shape_value(i, q) *
                                  fe_values.JxW(q);

              // Diffusion term
              cell_residual(i) -= spreading_coefficient *
                                  solution_gradient_loc[q] *
                                  fe_values.shape_grad(i, q) * fe_values.JxW(q);

              // Growth term
              cell_residual(i) += ((growth_coefficient * solution_loc[q]) *
                                   (1 - solution_loc[q])) *
                                  fe_values.shape_value(i, q) *
                                  fe_values.JxW(q);

              // Forcing term equals to 0, no contribution to cell_residual(i)
            }

          cell->get_dof_indices(dof_indices);

          jacobian_matrix.add(dof_indices, cell_matrix);
          residual_vector.add(dof_indices, cell_residual);
        }
    }
}

void
FisherKolmogorov1D::solve_linear_system()
{
  // setting for solver
  SolverControl solver_control(
    1000, 1e-6 * residual_vector.l2_norm()); // residual norm taken into account
                                             // for a more reliable tolerance

  // GMRES solver
  SolverGMRES<Vector<double>> solver(solver_control);

  std::cout << "  Solving the linear system" << std::endl;

  // solve with GMRES solver
  // TODO PRECONDITION
  solver.solve(jacobian_matrix, delta, residual_vector, PreconditionIdentity());
  std::cout << "  " << solver_control.last_step() << " GMRES iterations"
            << std::endl;
}

void
FisherKolmogorov1D::solve_newton()
{
  // parameters for the GMRES solver
  const unsigned int n_max_iters        = 1000;
  const double       residual_tolerance = 1e-6;
  unsigned int       n_iter             = 0;
  double             residual_norm      = residual_tolerance + 1;

  // Apply the boundary conditions: in this case no need for dirichlet boundary
  // conditions

  // iterative cycle till convergence
  while (n_iter < n_max_iters && residual_norm > residual_tolerance)
    {
      // assembly of the jacobian matrix and the residual vector
      assemble_system();
      residual_norm = residual_vector.l2_norm();

      std::cout << "  Newton iteration " << n_iter << "/" << n_max_iters
                << " - ||r|| = " << std::scientific << std::setprecision(6)
                << residual_norm << std::flush;

      // We actually solve the system only if the residual is larger than the
      // tolerance.
      if (residual_norm > residual_tolerance)
        {
          solve_linear_system();

          solution += delta;
        }
      else
        {
          std::cout << " < tolerance" << std::endl;
        }

      ++n_iter;
    }
}

void
FisherKolmogorov1D::solve()
{
  std::cout << "===============================================" << std::endl;

  time = 0.0;

  // Apply the initial condition.
  {
    std::cout << "Applying the initial condition" << std::endl;

    VectorTools::interpolate(dof_handler, c_0, solution);

    // Output the initial solution.
    output(0);
    std::cout << "-----------------------------------------------" << std::endl;
  }

  unsigned int time_step = 0;

  while (
    time <
    T -
      0.5 *
        deltat) // - 0.5*deltat is a tolerance setted to avoid extra-iterations
    {
      time += deltat;
      ++time_step;

      // Store the old solution, so that it is available for assembly.
      solution_old = solution;

      std::cout << "n = " << std::setw(3) << time_step
                << ", t = " << std::setw(5) << std::fixed << time << std::endl;

      // At every time step Newton's method to solve the non-linear system
      solve_newton();

      output(time_step);

      std::cout << std::endl;
    }
}


void
FisherKolmogorov1D::output(const unsigned int &time_step) const
{
  std::cout << "===============================================" << std::endl;
  std::cout << "Time step:  " << time_step << std::endl;

  // The DataOut class manages writing the results to a file.
  DataOut<dim> data_out;

  // It can write multiple variables (defined on the same mesh) to a single
  // file. Each of them can be added by calling add_data_vector, passing the
  // associated DoFHandler and a name.
  data_out.add_data_vector(dof_handler, solution, "solution");

  // Once all vectors have been inserted, call build_patches to finalize the
  // DataOut object, preparing it for writing to file.
  data_out.build_patches();

  // Then, use one of the many write_* methods to write the file in an
  // appropriate format.
  const std::string output_file_name =
    "output-" + std::to_string(time_step) + ".vtk";
  std::ofstream output_file(output_file_name);
  data_out.write_vtk(output_file);

  std::cout << "Output written to " << output_file_name << std::endl;

  std::cout << "===============================================" << std::endl;
}
