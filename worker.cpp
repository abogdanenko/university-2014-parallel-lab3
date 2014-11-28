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
    if (!args.usage_flag)
    {
        zero_matrix(U, n);
        zero_matrix(U_next, n);
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

void Worker::SetTopBC()
{
    Vector& row = U.back();
    for (Index j = 1; j < n - 1; j++)
    {
        const double x = Coord(j);
        row[j] = sin(M_PI * x) * exp(-1.0 * x);
    }
    U_next.back() = row;
}

void Worker::SetBottomBC()
{
    Vector& row = U.front();
    for (Index j = 1; j < n - 1; j++)
    {
        const double x = Coord(j);
        row[j] = sin(M_PI * x);
    }
    U_next.front() = row;
}

void Worker::MatrixWriteToFile() const
{
    ofstream fs;
    ostream& s = (args.matrix_filename == "-") ? cout :
        (fs.open(args.matrix_filename.c_str()), fs);

    s << fixed << setprecision(4);
    for (Index i = 0; i < n; i++)
    {
        for (Index j = 0; j < n; j++)
        {
            s << U_next[i][j] << ' ';
        }

        s << endl;
    }
}

void Worker::StatsWriteToFile() const
{
    ofstream fs;
    ostream& s = (args.stats_filename == "-") ? cout :
        (fs.open(args.stats_filename.c_str()), fs);

    const int w = string("iterations max: ").length();


    s
        << left
        << setw(w)
        << "matrix size:"
        << n
        << "x"
        << n
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
            const double diff = U[i][j] - U_next[i][j];
            const double abs_diff = fabs(diff);
            if (abs_diff > epsilon)
            {
                epsilon = abs_diff;
            }
        }
    }
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

void Worker::CalculateNextMatrix()
{
    if (args.omega == 0.0)
    {
        CalculateOmega();
    }

    for (Index i = 1; i < n - 1; i++)
    {
        for (Index j = 1; j < n - 1; j++)
        {
            const double p = U[i][j];
            const double q = 0.25 * (
                U[i - 1][j] +
                U[i + 1][j] +
                U[i][j - 1] +
                U[i][j + 1]);
            U_next[i][j] = omega * q + (1.0 - omega) * p;
        }
    }
}

void Worker::LogWrite() const
{
    if (args.verbose_flag)
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
    CalculateNextMatrix();
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
        Parser::PrintUsage();
    }
    else
    {
        timer.Start();
        SetTopBC();
        SetBottomBC();
        Loop();
        timer.Stop();
        MatrixWriteToFile();
        StatsWriteToFile();
    }
}
