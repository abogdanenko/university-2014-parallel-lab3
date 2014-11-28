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

void zero_matrix(Matrix& A, const Index size)
{
    return zero_matrix(A, size, size);
}

double GetWallTime()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return (double) time.tv_sec + (double) time.tv_usec * 1.0e-6;
}
