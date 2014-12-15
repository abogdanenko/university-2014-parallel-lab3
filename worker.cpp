#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <omp.h>

#include "worker.h"
#include "parser.h"
#include "routines.h"

using std::ofstream;
using std::ostream;

using std::cout;
using std::endl;
using std::setw;
using std::left;

using std::fixed;
using std::scientific;
using std::setprecision;

Worker::Worker(const Args& args):
    n(args.n),
    args(args),
    omega(1.0)
{
    #pragma omp parallel
    {
        #pragma omp single
        nt = omp_get_num_threads();
    }

    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    factor3(np, npx, npy, npz);

    nx = n / npx;
    ny = n / npy;
    nz = n / npz;

    vector<int> dims(3);
    vector<int> periods(3, 0);
    vector<int> coords(3);

    dims[0] = npx;
    dims[1] = npy;
    dims[2] = npz;

    MPI_Cart_create(MPI_COMM_WORLD, 3, &dims[0], &periods[0], 1, &comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Cart_coords(comm, rank, 3, &coords[0]);

    px = coords[0];
    py = coords[1];
    pz = coords[2];

    master_rank = Rank(0, 0, 0);
    is_master = rank == master_rank;

    if (!args.usage_flag)
    {
        zero_array3d(U, nx, ny, nz);
        zero_array3d(U_next, nx, ny, nz);
        H.Init(nx, ny, nz);
        if (args.omega != 0.0)
        {
            omega = args.omega;
        }
    }
}

double Worker::Coord(const Index i) const
{
    const double coef = 1.0 / (n - 1);
    return coef * i;
}

void Worker::SetBC()
{
    Matrix& side = U.back();
    for (Index j = 0; j < ny; j++)
    {
        const Index J = py * ny + j;
        if (J == 0 || J == n - 1)
        {
            continue;
        }
        const double y = Coord(J);

        for (Index k = 0; k < nz; k++)
        {
            const Index K = pz * nz + k;
            if (K == 0 || K == n - 1)
            {
                continue;
            }
            const double z = Coord(K);
            side[j][k] = border_condition(y, z);
        }
    }
    U_next.back() = side;
}

void Worker::ArrayWriteToFile() const
{
    if (args.array_filename.empty())
    {
        return;
    }

    if (is_master)
    {
        Vector buf(nz);

        ofstream fs;
        ostream& s = (args.array_filename == "-") ? cout :
            (fs.open(args.array_filename.c_str()), fs);

        s << fixed << setprecision(4);

        for (int pi = 0; pi < npx; pi++)
        {
            for (Index i = 0; i < nx; i++)
            {
                for (int pj = 0; pj < npy; pj++)
                {
                    for (Index j = 0; j < ny; j++)
                    {
                        for (int pk = 0; pk < npz; pk++)
                        {
                            MPI_Barrier(comm);
                            if (pi == 0 && pj == 0 && pk == 0)
                            {
                                buf = U_next[i][j];
                            }
                            else
                            {
                                const int worker = Rank(pi, pj, pk);
                                MPI_Recv(
                                    &buf[0],
                                    buf.size(),
                                    MPI_DOUBLE,
                                    worker,
                                    MPI_ANY_TAG,
                                    comm,
                                    MPI_STATUS_IGNORE);
                            }
                            vector_write_to_stream(s, buf);
                            MPI_Barrier(comm);
                        }
                        s << endl;
                    }
                }
            }
        }
    }
    else
    {
        for (int pi = 0; pi < npx; pi++)
        {
            for (Index i = 0; i < nx; i++)
            {
                for (int pj = 0; pj < npy; pj++)
                {
                    for (Index j = 0; j < ny; j++)
                    {
                        for (int pk = 0; pk < npz; pk++)
                        {
                            MPI_Barrier(comm);
                            if (pi == px && pj == py && pk == pz)
                            {
                                Vector buf = U_next[i][j];
                                MPI_Send(
                                    &buf[0],
                                    buf.size(),
                                    MPI_DOUBLE,
                                    master_rank,
                                    0,
                                    comm);
                            }
                            MPI_Barrier(comm);
                        }
                    }
                }
            }
        }
    }
}

void Worker::StatsWriteToFile() const
{
    if (!is_master)
    {
        return;
    }

    ofstream fs;
    ostream& s = (args.stats_filename == "-") ? cout :
        (fs.open(args.stats_filename.c_str()), fs);

    const int w = string("iterations max: ").length();

    s
        << left
        << setw(w)
        << "array size:"
        << n
        << "x"
        << n
        << "x"
        << n
        << endl;
    s
        << left
        << setw(w)
        << "array slice:"
        << nx
        << "x"
        << ny
        << "x"
        << nz
        << endl;
    s
        << left
        << setw(w)
        << "processes:"
        << np
        << endl;
    s
        << left
        << setw(w)
        << "topology:"
        << npx
        << "x"
        << npy
        << "x"
        << npz
        << endl;
    s
        << left
        << setw(w)
        << "threads:"
        << nt
        << endl;
    s
        << left
        << setw(w)
        << "epsilon:"
        << scientific
        << setprecision(3)
        << epsilon
        << endl;
    s
        << left
        << setw(w)
        << "epsilon max:"
        << args.epsilon_max
        << endl;
    s
        << left
        << setw(w)
        << "iterations:"
        << iteration_count
        << endl;
    s
        << left
        << setw(w)
        << "iterations max:";
    if (args.iteration_count_max)
    {
        s  << args.iteration_count_max;
    }
    else
    {
      s  << "no limit";
    }

    s
        << endl;
    s
        << left
        << setw(w)
        << "wall time, s:"
        << timer.GetDelta()
        << endl;
}

void Worker::CalculateEpsilon()
{
    epsilon = 0.0;

    /* workaround for old compiler that doesn't support
    #pragma omp for reduction(max:eps) */

    #pragma omp parallel
    {
        double eps = 0.0;
        for (Index i = 0; i < nx; i++)
        {
            for (Index j = 0; j < ny; j++)
            {
                for (Index k = 0; k < nz; k++)
                {
                    const double diff = U[i][j][k] - U_next[i][j][k];
                    const double abs_diff = fabs(diff);
                    if (abs_diff > eps)
                    {
                        eps = abs_diff;
                    }
                }
            }
        }

        #pragma omp critical
        if (eps > epsilon)
        {
            epsilon = eps;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &epsilon, 1, MPI_DOUBLE, MPI_MAX, comm);
}

void Worker::CalculateOmega()
{
    const double c1 = M_PI_2 / (n + 1);
    const double rho = 1.0 - c1 * c1;
    const double rho2 = rho * rho;

    if (iteration_count == 1)
    {
        omega = 1.0 / (1.0 - rho2 / 2.0);
    }

    if (iteration_count > 1)
    {
        omega = 1.0 / (1.0 - rho2  * omega / 4.0);
    }
}

int Worker::Rank(const int pi, const int pj, const int pk) const
{
    int rank;
    vector<int> coords(3);
    coords[0] = pi;
    coords[1] = pj;
    coords[2] = pk;
    MPI_Cart_rank(comm, &coords[0], &rank);
    return rank;
}

void Worker::SendReceiveMatrix(Matrix& A, const int dest)
{
    for (Index i = 0; i < A.size(); i++)
    {
        Vector& buf = A[i];
        MPI_Sendrecv_replace(
            &buf[0],
            buf.size(),
            MPI_DOUBLE,
            dest,
            0,
            dest,
            MPI_ANY_TAG,
            comm,
            MPI_STATUS_IGNORE);
    }
}

void Worker::SendReceiveX0()
{
    if (px == 0)
    {
        return;
    }

    H.x0_exists = true;
    H.x0 = U.front();
    const int dest = Rank(px - 1, py, pz);
    SendReceiveMatrix(H.x0, dest);
}

void Worker::SendReceiveX1()
{
    if (px == npx - 1)
    {
        return;
    }

    H.x1_exists = true;
    H.x1 = U.back();
    const int dest = Rank(px + 1, py, pz);
    SendReceiveMatrix(H.x1, dest);
}

void Worker::SendReceiveY0()
{
    if (py == 0)
    {
        return;
    }

    H.y0_exists = true;
    for (Index i = 0; i < nx; i++)
    {
        H.y0[i] = U[i][0];
    }

    const int dest = Rank(px, py - 1, pz);
    SendReceiveMatrix(H.y0, dest);
}

void Worker::SendReceiveY1()
{
    if (py == npy - 1)
    {
        return;
    }

    H.y1_exists = true;
    for (Index i = 0; i < nx; i++)
    {
        H.y1[i] = U[i][ny - 1];
    }

    const int dest = Rank(px, py + 1, pz);
    SendReceiveMatrix(H.y1, dest);
}

void Worker::SendReceiveZ0()
{
    if (pz == 0)
    {
        return;
    }

    H.z0_exists = true;
    for (Index i = 0; i < nx; i++)
    {
        for (Index j = 0; j < ny; j++)
        {
            H.z0[i][j] = U[i][j][0];
        }
    }

    const int dest = Rank(px, py, pz - 1);
    SendReceiveMatrix(H.z0, dest);
}

void Worker::SendReceiveZ1()
{
    if (pz == npz - 1)
    {
        return;
    }

    H.z1_exists = true;
    for (Index i = 0; i < nx; i++)
    {
        for (Index j = 0; j < ny; j++)
        {
            H.z1[i][j] = U[i][j][nz - 1];
        }
    }

    const int dest = Rank(px, py, pz + 1);
    SendReceiveMatrix(H.z1, dest);
}

void Worker::RunInOrder(
    WorkerMemberFunction f1,
    WorkerMemberFunction f2,
    const bool order)
{
    if (order)
    {
        (this->*f1)();
        (this->*f2)();
    }
    else
    {
        (this->*f2)();
        (this->*f1)();
    }
}

void Worker::SendReceiveHalo()
{
    H.x0_exists = false;
    H.x1_exists = false;
    H.y0_exists = false;
    H.y1_exists = false;
    H.z0_exists = false;
    H.z1_exists = false;
    RunInOrder(&Worker::SendReceiveX0, &Worker::SendReceiveX1, px % 2);
    RunInOrder(&Worker::SendReceiveY0, &Worker::SendReceiveY1, py % 2);
    RunInOrder(&Worker::SendReceiveZ0, &Worker::SendReceiveZ1, pz % 2);
}

void Worker::CalculateUNext()
{
    if (args.omega == 0.0)
    {
        CalculateOmega();
    }

    #pragma omp parallel for
    for (Index i = 0; i < nx; i++)
    {
        if (i == 0 && !H.x0_exists)
        {
            continue;
        }

        if (i == nx - 1 && !H.x1_exists)
        {
            continue;
        }

        const Matrix& x0 = i == 0      ? H.x0 : U[i - 1];
        const Matrix& x1 = i == nx - 1 ? H.x1 : U[i + 1];

        HaloSection h;
        h.TakeSection(H, i);
        CalculateSection(U_next[i], U[i], x0, x1, h);
    }
}

void Worker::CalculateSection(
    Matrix& next,
    const Matrix& current,
    const Matrix& mx0,
    const Matrix& mx1,
    const HaloSection& h)
{
    for (Index j = 0; j < ny; j++)
    {
        if (j == 0 && !h.y0_exists)
        {
            continue;
        }

        if (j == ny - 1 && !h.y1_exists)
        {
            continue;
        }

        for (Index k = 0; k < nz; k++)
        {
            if (k == 0 && !h.z0_exists)
            {
                continue;
            }

            if (k == nz - 1 && !h.z1_exists)
            {
                continue;
            }

            const double x0  = mx0[j][k];
            const double x1  = mx1[j][k];
            const double y0  = j == 0      ? h.y0[k] : current[j - 1][k];
            const double y1  = j == ny - 1 ? h.y1[k] : current[j + 1][k];
            const double z0  = k == 0      ? h.z0[j] : current[j][k - 1];
            const double z1  = k == nz - 1 ? h.z1[j] : current[j][k + 1];
            const double sum = x0 + x1 + y0 + y1 + z0 + z1;
            const double p   = current[j][k];
            const double q   = sum / 6.0;
            next[j][k]       = omega * q + (1.0 - omega) * p;
        }
    }
}

void Worker::LogWrite() const
{
    if (args.verbose_flag && is_master)
    {
        cout << "iteration:"
            << setw(8)
            << iteration_count
            << "    epsilon: "
            << scientific
            << setprecision(3)
            << epsilon
            << "    omega: "
            << setprecision(6)
            << omega
            << endl;
    }
}

void Worker::Step()
{
    U.swap(U_next);
    SendReceiveHalo();
    CalculateUNext();
    CalculateEpsilon();
    LogWrite();
    iteration_count++;
}

void Worker::Loop()
{
    iteration_count = 0;

    if (args.iteration_count_max)
    {
        do
        {
            Step();
        }
        while (epsilon > args.epsilon_max &&
            iteration_count < args.iteration_count_max);
    }
    else
    {
        do
        {
            Step();
        }
        while (epsilon > args.epsilon_max);
    }

}

void Worker::Run()
{
    if (args.usage_flag)
    {
        if (is_master)
        {
            Parser::PrintUsage();
        }
    }
    else
    {
        MPI_Barrier(comm);
        timer.Start();
        if (px == npx - 1)
        {
            SetBC();
        }
        Loop();
        MPI_Barrier(comm);
        timer.Stop();
        ArrayWriteToFile();
        StatsWriteToFile();
    }
}
