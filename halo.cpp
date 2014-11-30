#include "halo.h"
#include "routines.h"

void Halo::Init(const int nx, const int ny, const int nz)
{
    zero_matrix(x0, ny, nz);
    zero_matrix(x1, ny, nz);
    zero_matrix(y0, nx, nz);
    zero_matrix(y1, nx, nz);
    zero_matrix(z0, nx, ny);
    zero_matrix(z1, nx, ny);
}
