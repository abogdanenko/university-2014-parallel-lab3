#ifndef WORKER_H
#define WORKER_H

#include "args.h"
#include "typedefs.h"
#include "timer.h"

class Worker;
typedef void (Worker::* WorkerMemberFunction)();

class Worker
{
    // mesh size
    Index n;
    Index nx;
    Index ny;
    Index nz;
    // coords of this processor
    int px;
    int py;
    int pz;
    // total number of processors and number of processors in each dimension
    int np;
    int npx;
    int npy;
    int npz;
    MPI_Comm comm;
    bool is_master;
    int master_rank;
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
    int Rank(const int pi, const int pj, const int pk) const;
    void SendReceiveMatrix(Matrix& A, const int dest);
    void SendReceiveX0();
    void SendReceiveX1();
    void SendReceiveY0();
    void SendReceiveY1();
    void SendReceiveZ0();
    void SendReceiveZ1();
    void SendReceiveHalo();
    void RunInOrder(
        WorkerMemberFunction f1,
        WorkerMemberFunction f2,
        const bool order);
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
