#include "halosection.h"

void HaloSection::TakeSection(const Halo& H, const int i)
{
    y0_exists = H.y0_exists;
    y1_exists = H.y1_exists;
    z0_exists = H.z0_exists;
    z1_exists = H.z1_exists;

    if (y0_exists)
    {
        y0 = H.y0[i];
    }

    if (y1_exists)
    {
        y1 = H.y1[i];
    }

    if (z0_exists)
    {
        z0 = H.z0[i];
    }

    if (z1_exists)
    {
        z1 = H.z1[i];
    }
}
