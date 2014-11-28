#include <sys/time.h>

#include "routines.h"

double triangle(const double x)
{
    return 1.0 - abs(2.0 * x - 1.0);
}

double border_condition(const double y, const double z)
{
    return triangle(y) * triangle(z);
}

void zero_matrix(Matrix& A, const Index rows, const Index cols)
{
    A.resize(rows);
    for (Index i = 0; i < rows; i++)
    {
        Vector& row = A[i];
        row.resize(cols);
        fill(row.begin(), row.end(), 0.0);
    }
}

void zero_array3d(Array3d& A, const Index n1, const Index n2, const Index n3)
{
    A.resize(n1);
    for (Index i = 0; i < n1; i++)
    {
        Matrix& M = A[i];
        zero_matrix(M, n2, n3);
    }
}

double GetWallTime()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return (double) time.tv_sec + (double) time.tv_usec * 1.0e-6;
}
