#include "FisherKolmogorov.hpp"

void 
FisherKolmogorov::setup() 
{
  Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);
}

void
HeatNonLinear::solve_linear_system()
{
  //setting for solver
  SolverControl solver_control(1000, 1e-6 * residual_vector.l2_norm()); //residual norm taken into account for a more reliable tolerance

  //GMRES solver 
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);

  //AMG preconditioner declaration and initialization
  TrilinosWrappers::PreconditionAMG preconditioner;
  preconditioner.initialize(jacobian_matrix, TrilinosWrappers::PreconditionAMG::AdditionalData(1.0));
  
  //solve with GMRES solver
  solver.solve(jacobian_matrix, delta_owned, residual_vector, preconditioner);
  pcout << "  " << solver_control.last_step() << " CG iterations" << std::endl;
}

void
FisherKolmogorov::solve_newton()
{
  //parameters for the GMRES solver
  const unsigned int n_max_iters        = 1000;
  const double       residual_tolerance = 1e-6;
  unsigned int n_iter        = 0;
  double       residual_norm = residual_tolerance + 1;

  //Apply the boundary conditions: in this case no need for dirichlet boundary conditions


  //iterative cycle till convergence
  while (n_iter < n_max_iters && residual_norm > residual_tolerance)
    {
      //assembly of the jacobian matrix and the residual vector
      assemble_system();
      residual_norm = residual_vector.l2_norm();

      pcout << "  Newton iteration " << n_iter << "/" << n_max_iters
            << " - ||r|| = " << std::scientific << std::setprecision(6)
            << residual_norm << std::flush;

      // We actually solve the system only if the residual is larger than the tolerance.
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

  while (time < T - 0.5 * deltat) // - 0.5*deltat is a tolerance setted to avoid extra-iterations
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


//output function
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
