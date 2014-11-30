#ifndef HALOSECTION_H
#define HALOSECTION_H

#include "halo.h"

class HaloSection
{
    public:
    bool y0_exists;
    bool y1_exists;
    bool z0_exists;
    bool z1_exists;

    Vector y0;
    Vector y1;
    Vector z0;
    Vector z1;

    void TakeSection(const Halo& H, const int i);
};

#endif
