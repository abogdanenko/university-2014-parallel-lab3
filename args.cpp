#include "args.h"

Args::Args():
    n(512),
    epsilon_max(0.01),
    omega(1.0),
    iteration_count_max(0),
    usage_flag(false),
    verbose_flag(false),
    matrix_filename("-"),
    stats_filename("-")
{

}
