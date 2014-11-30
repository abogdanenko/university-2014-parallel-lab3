#include <sys/time.h>
#include <cmath>

#include "routines.h"

double triangle(const double x)
{
    return 1.0 - fabs(2.0 * x - 1.0);
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

void vector_write_to_stream(ostream& s, const Vector& v)
{
    for (Index i = 0; i < v.size(); i++)
    {
        s << v[i] << ' ';
    }
}

void factor3(int n, int& x, int& y, int& z)
{
    x = 1;
    y = 1;
    z = 1;

    const vector<int*> p = {&x, &y, &z};
    int i = 0;
    while (n > 1)
    {
        n /= 2;
        *p[i] *= 2;
        i = (i + 1) % 3;
    }
}
