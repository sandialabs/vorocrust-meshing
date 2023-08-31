#include "Config.h"
#include <iostream>

/*
 * This test is only enabled if VOROCRUST_ENABLE_MPI is set, so if we're
 * building this test, we better have USE_MPI defined. To test this we
 * just create two mains. One will exit with an error if MPI doesn't exist
 * and the other is a standard "Hello World" style MPI app.
 */
#if defined USE_MPI


#include <mpi.h>

int main(int argc, char* argv[])
{
    int rank = -1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::cout << "Hello from node rank " << rank << std::endl;

    MPI_Finalize();
    return 0;
}


#else


int main(int argc, char* argv[])
{
    std::cout << "ERROR! USE_MPI was not defined, but this test requires MPI" << std::endl;
    return 1;
}


#endif
