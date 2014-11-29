#include "worker.h"
#include "parser.h"
#include "routines.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>

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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    nx = n / px;
    ny = n / py;
    nz = n / pz;

    if (!args.usage_flag)
    {
        zero_array3d(U, nx, ny, nz);
        zero_array3d(U_next, nx, ny, nz);
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
    for (Index j = 1; j < n - 1; j++)
    {
        for (Index k = 1; k < n - 1; k++)
        {
            const double y = Coord(j);
            const double z = Coord(k);
            side[j][k] = border_condition(y, z);
        }
    }
    U_next.back() = side;
}

void Worker::ArrayWriteToFile() const
{
    if (master)
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
                                int worker;
                                const vector<int> coords = {pi, pj, pk};
                                MPI_Cart_rank(comm, &coords[0], &worker);
                                MPI_recv(
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
                                const Vector& buf = U_next[i][j];
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
    if (rank)
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
        << world_size
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
    for (Index i = 0; i < n; i++)
    {
        for (Index j = 0; j < n; j++)
        {
            for (Index k = 0; k < n; k++)
            {
                const double diff = U[i][j][k] - U_next[i][j][k];
                const double abs_diff = fabs(diff);
                if (abs_diff > epsilon)
                {
                    epsilon = abs_diff;
                }
            }
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

void Worker::CalculateUNext()
{
    if (args.omega == 0.0)
    {
        CalculateOmega();
    }

    for (Index i = 0; i < nx; i++)
    {
        if (i == 0 && !O3d.front_exists || i == nx - 1 && !O3d.back_exists)
        {
            continue;
        }

        const Matrix& front = i == 0 ? O3d.front : U_next[i - 1];
        const Matrix& back i == nx - 1 ? O3d.back : U_next[i + 1];

        Neighborhood O(front, back);

        if (O3d.left_exists)
        {
            O.SetLeft(O3d.left[i]);
        }

        if (O3d.right_exists)
        {
            O.SetRight(O3d.right[i]);
        }

        if (O3d.bottom_exists)
        {
            O.SetBottom(O3d.bottom[i]);
        }

        if (O3d.top_exists)
        {
            O.SetTop(O3d.top[i]);
        }

        CalculateMatrix(U_next[i], U[i], O);
    }
}

void Worker::CalculateMatrix(
    Matrix& next,
    const Matrix& current,
    const Neighborhood& O)
{
    for (Index j = 0; j < ny; j++)
    {
        if (j == 0 && !O.left_exists || j == ny - 1 && !O.right_exists)
        {
            continue;
        }

        for (Index k = 0; k < nz; k++)
        {
            if (k == 0 && !O.bottom_exists || k == nz - 1 && !O.top_exists)
            {
                continue;
            }

            const double left   = j == 0      ? O.left[k]   : current[j - 1][k];
            const double right  = j == ny - 1 ? O.right[k]  : current[j + 1][k];
            const double bottom = k == 0      ? O.bottom[j] : current[j][k - 1];
            const double top    = k == nz - 1 ? O.top[j]    : current[j][k + 1];
            const double back   = O.back[j][k];
            const double front  = O.front[j][k];
            const double sum    = front + back + left + right + bottom + top;
            const double p      = current[j][k];
            const double q      = sum / 6.0;
            next[j][k]          = omega * q + (1.0 - omega) * p;
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
        if (rank == 0)
        {
            Parser::PrintUsage();
        }
    }
    else
    {
        MPI_Barrier(MPI_COMM_WORLD);
        timer.Start();
        if (pi == px - 1)
        {
            SetBC();
        }
        Loop();
        MPI_Barrier(MPI_COMM_WORLD);
        timer.Stop();
        ArrayWriteToFile();
        StatsWriteToFile();
    }
}
