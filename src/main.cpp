#include "FisherKolmogorov.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

int main (int argc, char *argv[]) {
    
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    const std::string  mesh_filename = "../mesh/brain-h3.0.msh";
    const unsigned int degree        = 1;

    FisherKolmogorov problem(mesh_filename, degree);
    
}

