#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "FisherKolmogorov.hpp"

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string  mesh_filename = "../brain-h3.0.msh";
  const unsigned int degree        = 1;
  const double       deltat        = 0.1;
  const double       T             = 1;
  const double       dext          = 0.0001;

  FisherKolmogorov problem(mesh_filename, degree, deltat, T, dext);
  problem.setup();
  problem.solve();
}
