#ifndef WORKER_H
#define WORKER_H

#include "args.h"
#include "typedefs.h"
#include "timer.h"

class Worker
{
    // mesh size
    Index n;
    Index nx;
    Index ny;
    Index nz;
    // coords of this processor
    int rank;
    int px;
    int py;
    int pz;
    // total number of processors and number of processors in each dimension
    int world_size;
    int npx;
    int npy;
    int npz;
    Args args;
    Timer timer;
    Array3d U;
    Array3d U_next;
    Halo H;
    double epsilon;
    unsigned iteration_count;
    double omega;

    double Coord(const Index i) const;
    void SetBC();
    void LogWrite() const;
    void ArrayWriteToFile() const;
    void StatsWriteToFile() const;
    void CalculateEpsilon();
    void CalculateOmega();
    void CalculateSection(
        Matrix& next,
        const Matrix& current,
        const Matrix& mx0,
        const Matrix& mx1,
        const HaloSection& h);
    void CalculateUNext();
    void Step();
    void Loop();
    public:
    void Run();
    Worker(const Args& args);
};

#endif
