# Parallel Jacobi Laplace 3D Solver
Hybrid (MPI+Openmp) parallel program that solves 3D laplace equation using
Jacobi method.

Number of points, iteration limit, output file name etc. should be supplied
as command line arguments. See *Parser::PrintUsage* function in file
*parser.cpp* for program usage.

The cube-shaped domain is divided into n1 x n2 x n3 = np smaller regions which
are distributed between workers. Workers exchange border regions using messages
(MPI).

Each worker process solves the equation in its region in parallel using threads
(OpenMP).

When required tolerance or maximum number of iterations is reached the
distributed array is written to file by master process.

Number of iterations, time and other computation details are also recorded.
