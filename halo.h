#ifndef HALO_H
#define HALO_H

#include "typedefs.h"

class Halo
{
    public:
    bool x0_exists;
    bool x1_exists;
    bool y0_exists;
    bool y1_exists;
    bool z0_exists;
    bool z1_exists;

    Matrix x0;
    Matrix x1;
    Matrix y0;
    Matrix y1;
    Matrix z0;
    Matrix z1;

    void Init(const int nx, const int ny, const int nz);
};

#endif
