#ifndef ARGS_H
#define ARGS_H

#include <string>

using std::string;

class Args
{
    public:
    unsigned n;
    double epsilon_max;
    double omega;
    unsigned iteration_count_max;
    bool usage_flag;
    bool verbose_flag;
    // "-" means 'write to stdout' (default)
    string matrix_filename;
    string stats_filename;

    Args();
};

#endif
