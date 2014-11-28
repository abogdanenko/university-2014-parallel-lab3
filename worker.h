#ifndef WORKER_H
#define WORKER_H

#include "args.h"
#include "typedefs.h"
#include "timer.h"

class Worker
{
    Index n;
    Args args;
    Timer timer;
    Array3d U;
    Array3d U_next;
    double epsilon;
    unsigned iteration_count;
    double omega;
    double Coord(const Index i) const;
    void SetBC();
    void LogWrite() const;
    void MatrixWriteToFile() const;
    void StatsWriteToFile() const;
    void CalculateEpsilon();
    void CalculateOmega();
    void CalculateNextMatrix();
    void Step();
    void Loop();
    public:
    void Run();
    Worker(const Args& args);
};

#endif
