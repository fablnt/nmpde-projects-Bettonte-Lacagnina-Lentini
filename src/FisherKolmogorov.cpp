#include "../include/FisherKolmogorov.hpp"

/**
 * Function called by main() function.
 * Creates the setup of the problem.
 * Takes the mesh as input, requesting an acceptable .msh file from deal.II;
 * in case the mesh is built from scratch (as in our 2D case) it is necessary to
 * follow the 'step-49 tutorial' in the reference documentation for deal.II
 * (referring to its version 9.5.0). Initializes the finite element space
 * choosing the 'Fe_SimplexP' in order to deal with 2D or 3D mesh with
 * tetrahedral elements and chooses the Gauss-Legendre quadrature formula to
 * compute the required integrals. Initializes the Dof Handler, to identify and
 * interact with the degrees of freedom located in the mesh. Initializes the
 * grey and white matter region in order to set different parameters depending
 * on the type of element considered.
 * Computes the center of the domain, necessary to implement the axonal
 * transport coefficient (parameter), and initializes the fiber orientation .
 * Initializes the structure of the linear system, which is what the starting
 * problem is reduced to.
 */
void
FisherKolmogorov::setup()
{
  Triangulation<dim> mesh_serial;

  GridIn<dim> grid_in;
  grid_in.attach_triangulation(mesh_serial);

  std::ifstream grid_in_file(mesh_file_name);
  grid_in.read_msh(grid_in_file);

  GridTools::partition_triangulation(mpi_size, mesh_serial);
  const auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      mesh_serial, MPI_COMM_WORLD);
  mesh.create_triangulation(construction_data);


  pcout << "  Number of elements = " << mesh.n_global_active_cells()
        << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_SimplexP<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;
  }

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  // Initialize grey and white matter regions by setting the material_id to 1
  // for grey matter elements and leaving it to 0 for white matter ones.
  {
    pcout << "Initializing the white and grey matter regions" << std::endl;
    grey_matter = Grey_matter<dim>();

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;

        if (grey_matter.check_grey_matter(cell->center()))
          cell->set_material_id(1);
      }
  }

  // Compute the center of the domain.
  {
    Point<dim> local_center;
    local_center.clear();
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;
        local_center += cell->center();
      }

    MPI_Allreduce(
      &local_center, &global_center, dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    global_center /= mesh.n_global_active_cells();

    pcout << "Center of the domain = " << global_center << std::endl;
  }

  // Initialize the fiber orientation.
  {
    pcout << "Initializing the fiber orientation" << std::endl;
    pcout << "  Orientation = " << orientation << std::endl;
    direction = get_direction<dim>(orientation, global_center);
  }

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);
    sparsity.compress();

    pcout << "  Initializing the matrices" << std::endl;
    jacobian_matrix.reinit(sparsity);

    pcout << "  Initializing the system right-hand side" << std::endl;
    residual_vector.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    delta_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    solution_old = solution;
  }
}

/**
 * Assembles the linear system associated to the tangent problem (Newton
 * method). For all elements (cells) of the domain (mesh) does two operation:
 * setting of the specific parameters and calculation of their local
 * contribution in terms of linear system to the final system that as to be
 * solved. The setting of the parameters is done in relation to two aspects:
 * type of diffusion (isotropic, circumferential, radial and axon-based fiber
 * orientations) and the type of element (grey or white matter). When this
 * process is finished, all active MPI processes combine their local system to
 * obtain the final one to solve.
 */
void
FisherKolmogorov::assemble_system()
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

  double         growth_coefficient_loc;
  Tensor<2, dim> spreading_coefficient_loc;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_matrix   = 0.0;
      cell_residual = 0.0;

      fe_values.get_function_values(solution, solution_loc);
      fe_values.get_function_gradients(solution, solution_gradient_loc);
      fe_values.get_function_values(solution_old, solution_old_loc);

      for (unsigned int q = 0; q < n_q; ++q)
        {
          Vector<double> n_loc(dim);
          direction->vector_value(fe_values.quadrature_point(q), n_loc);

          // Convert the direction vector to a Tensor.
          Tensor<1, dim> n_loc_tensor;
          for (unsigned int i = 0; i < dim; ++i)
            n_loc_tensor[i] = n_loc[i];

          // Set the growth and spreading coefficients depending on the type
          // of element: grey matter (material_id = 1) or white matter
          // (material_id = 0).
          if (cell->material_id() == 1)
            {
              growth_coefficient_loc =
                growth_coefficient_grey.value(fe_values.quadrature_point(q));

              spreading_coefficient_loc =
                isotropic_coefficient_grey.value(
                  fe_values.quadrature_point(q)) *
                  unit_symmetric_tensor<dim>() +
                anisotropic_coefficient_grey.value(
                  fe_values.quadrature_point(q)) *
                  outer_product(n_loc_tensor, n_loc_tensor);
            }

          else
            {
              growth_coefficient_loc =
                growth_coefficient_white.value(fe_values.quadrature_point(q));

              spreading_coefficient_loc =
                isotropic_coefficient_white.value(
                  fe_values.quadrature_point(q)) *
                  unit_symmetric_tensor<dim>() +
                anisotropic_coefficient_white.value(
                  fe_values.quadrature_point(q)) *
                  outer_product(n_loc_tensor, n_loc_tensor);
            }

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
                    scalar_product(spreading_coefficient_loc *
                                     fe_values.shape_grad(i, q),
                                   fe_values.shape_grad(j, q)) *
                    fe_values.JxW(q);

                  // Non-linear stiffness matrix, second term.
                  cell_matrix(i, j) -=
                    growth_coefficient_loc * fe_values.shape_value(i, q) *
                    fe_values.shape_value(j, q) * fe_values.JxW(q);

                  // Non-linear stiffness matrix, third term.
                  cell_matrix(i, j) +=
                    2 * growth_coefficient_loc * solution_loc[q] *
                    fe_values.shape_value(i, q) * fe_values.shape_value(j, q) *
                    fe_values.JxW(q);
                }

              // Assemble the residual vector (with changed sign).
              // Time derivative term.
              cell_residual(i) -= (solution_loc[q] - solution_old_loc[q]) /
                                  deltat * fe_values.shape_value(i, q) *
                                  fe_values.JxW(q);

              // Diffusion term
              cell_residual(i) -= scalar_product(spreading_coefficient_loc *
                                                   solution_gradient_loc[q],
                                                 fe_values.shape_grad(i, q)) *
                                  fe_values.JxW(q);

              // Growth term
              cell_residual(i) += ((growth_coefficient_loc * solution_loc[q]) *
                                   (1 - solution_loc[q])) *
                                  fe_values.shape_value(i, q) *
                                  fe_values.JxW(q);


              // Forcing term equals to 0, no contribution to
              // cell_residual(i)
            }
        }
      cell->get_dof_indices(dof_indices);

      jacobian_matrix.add(dof_indices, cell_matrix);
      residual_vector.add(dof_indices, cell_residual);
    }
  // Parallel communication among processes
  jacobian_matrix.compress(VectorOperation::add);
  residual_vector.compress(VectorOperation::add);
}

/**
 * Solves the linear system obtained in the assemble_system() function,
 * parallelizing the computation on the different active MPI processes.
 * Uses the GMRES method as iterative method, setting the maximum number of
 * iterations equals to 1000, the tolerance equals to 1.e-6 and the Algebraic
 * Multigrid as preconditioner.
 */
void
FisherKolmogorov::solve_linear_system()
{
  // setting for solver
  SolverControl solver_control(
    1000,
    1e-6 * residual_vector.l2_norm()); // residual norm taken into account
                                       // for a more reliable tolerance

  // GMRES solver
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);

  // AMG preconditioner declaration and initialization
  TrilinosWrappers::PreconditionAMG preconditioner;
  preconditioner.initialize(
    jacobian_matrix, TrilinosWrappers::PreconditionAMG::AdditionalData(1.0));

  // solve with GMRES solver
  solver.solve(jacobian_matrix, delta_owned, residual_vector, preconditioner);

  pcout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;
}

/**
 * Solves the non-linear problem using the Newton method.
 * Invoked for each time step by the solve() function, finds a numerical
 * solution of the non-linear system assembling and solving a linear system,
 * using as maxinum number of iterations and as tolerance the same used for
 * the GMRES solver.
 */
void
FisherKolmogorov::solve_newton()
{
  // parameters for the GMRES solver
  const unsigned int n_max_iters        = 1000;
  const double       residual_tolerance = 1e-6;
  unsigned int       n_iter             = 0;
  double             residual_norm      = residual_tolerance + 1;

  // There are no dirichlet boundary conditions


  // iterative cycle till convergence
  while (n_iter < n_max_iters && residual_norm > residual_tolerance)
    {
      // assembly of the jacobian matrix and the residual vector
      assemble_system();
      residual_norm = residual_vector.l2_norm();

      pcout << "  Newton iteration " << n_iter << "/" << n_max_iters
            << " - ||r|| = " << std::scientific << std::setprecision(6)
            << residual_norm << std::flush;

      // We actually solve the system only if the residual is larger than the
      // tolerance.
      if (residual_norm > residual_tolerance)
        {
          solve_linear_system();

          solution_owned += delta_owned;
          solution = solution_owned;
        }
      else
        {
          pcout << " < tolerance" << std::endl;
        }

      ++n_iter;
    }
}

/**
 * Function called by main() function.
 * For each time step applies the Newton method (solving the tangent problem)
 * and calls the output() function which creates output files.
 */
void
FisherKolmogorov::solve()
{
  pcout << "===============================================" << std::endl;

  time = 0.0;

  // Apply the initial condition.
  {
    pcout << "Applying the initial condition" << std::endl;

    VectorTools::interpolate(dof_handler, c_0, solution_owned);
    solution = solution_owned;

    // Output the initial solution.
    output(0);
    pcout << "-----------------------------------------------" << std::endl;
  }

  unsigned int time_step = 0;

  while (time < T - 0.5 * deltat) // - 0.5*deltat is a tolerance setted to
                                  // avoid extra-iterations
    {
      time += deltat;
      ++time_step;

      // Store the old solution, so that it is available for assembly.
      solution_old = solution;

      pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
            << std::fixed << time << std::endl;

      // At every time step Newton's method to solve the non-linear system
      solve_newton();

      output(time_step);

      pcout << std::endl;
    }
}

/**
 * Output function.
 * Creates the output data for the selected time step, in particular a VTU
 * file for each MPI process, a PVTU file that points to VTU files and
 * contains information about splitting data across multiple MPI processes.
 * The output files are saved in the build directory.
 *
 * @param time_step   represents a distinct instant within the time interval T
 */
void
FisherKolmogorov::output(const unsigned int &time_step) const
{
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, solution, "u");

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  data_out.write_vtu_with_pvtu_record(
    "./", "output", time_step, MPI_COMM_WORLD, 3);
}
