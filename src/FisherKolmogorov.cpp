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
  SolverControl solver_control(1000, 1e-6 * residual_vector.l2_norm());

  ////what is the best solver methods? (to do)
  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

  //what is the best preconditioner?
  //I noticed that for other problems (lab8) AMG preconditioner works betetr that SSOR preconditioner
  //The actual reason why is still to determine properly 
  TrilinosWrappers::PreconditionAMG preconditioner;
  preconditioner.initialize(jacobian_matrix, TrilinosWrappers::PreconditionAMG::AdditionalData(1.0));
  solver.solve(jacobian_matrix, delta_owned, residual_vector, preconditioner);
  pcout << "  " << solver_control.last_step() << " CG iterations" << std::endl;
}

void
FisherKolmogorov::solve_newton()
{
  //parameters for the CG solver
  const unsigned int n_max_iters        = 1000;
  const double       residual_tolerance = 1e-6;
  unsigned int n_iter        = 0;
  double       residual_norm = residual_tolerance + 1;

  //Apply the boundary conditions (to do)
  //I am not yet confident on this step since the particular type of the mesh (brain)
  //and doubts about the dimension (dim = 1 or dim = 3)
  {
  }

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
          //why don't we compress here ? Does each process solve the complete system ?
          //propably since the solving step is managed by Trillinos, the compress is alrready done
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

  // Apply the initial condition (to do)
  //How does this step change according to different dimensions (1 or 3); still to figure out properly
  {
  }

  unsigned int time_step = 0;

  while (time < T - 0.5 * deltat)
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

