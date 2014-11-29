#include <cstdlib>
#include <iostream>
#include <mpi.h>

#include "parser.h"
#include "worker.h"

using std::cerr;
using std::endl;

int main(int argc, char** argv)
{
    int exit_code = EXIT_SUCCESS;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    try
    {
        Parser parser(argc, argv);
        Args args = parser.Parse();
        Worker worker(args);
        worker.Run();
    }
    catch (Parser::ParseError& e)
    {
        if (rank == 0)
        {
            cerr << e.what() << endl;
            Parser::PrintUsage();
        }
        exit_code = EXIT_FAILURE;
    }

    MPI_Finalize();
    return exit_code;
}
