#ifndef ROUTINES_H
#define ROUTINES_H

#include <sstream>

#include "typedefs.h"

using std::stringstream;
using std::string;

template <class Number>
Number string_to_number(const string& s)
{
    Number n;
    stringstream ss(s);
    ss >> n;
    return n;
}

double border_condition(const double y, const double z);
void zero_matrix(Matrix& A, const Index rows, const Index cols);
void zero_matrix(Matrix& A, const Index size);
double GetWallTime();

#endif
