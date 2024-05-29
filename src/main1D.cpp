#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "FisherKolmogorov1D.hpp"

int
main()
{
  const std::string  mesh_filename         = "../mesh1D.msh";
  const unsigned int degree                = 1;
  const double       deltat                = 0.1;
  const double       T                     = 20;
  const double       spreading_coefficient = 0.0002;
  const double       growth_coefficient    = 2;

  FisherKolmogorov1D problem(mesh_filename,
                             degree,
                             deltat,
                             T,
                             spreading_coefficient,
                             growth_coefficient);
  problem.setup();
  problem.solve();
}
