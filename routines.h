#ifndef ROUTINES_H
#define ROUTINES_H

#include <sstream>

#include "typedefs.h"

using std::stringstream;
using std::string;
using std::ostream;

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
void zero_array3d(Array3d& A, const Index n1, const Index n2, const Index n3);
void vector_write_to_stream(ostream& s, const Vector& v);
void factor3(const int n, int& x, int& y, int& z);

#endif
