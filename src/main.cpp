#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "FisherKolmogorov.hpp"

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string  mesh_filename  = "../mesh/brain-h3.0.msh";
  const unsigned int degree         = 1;
  const double       deltat         = 0.4;
  const double       T              = 40;
  const std::string  seeding_region = "Amyloid-Beta deposits";
  const std::string  orientation    = "axon-based";

  FisherKolmogorov problem(
    mesh_filename, degree, deltat, T, seeding_region, orientation);

  problem.setup();
  problem.solve();
}
