#ifndef PARSER_H
#define PARSER_H

#include <stdexcept>
#include <string>

#include "args.h"

using std::string;
using std::runtime_error;

class Parser
{
    const int argc;
    char** const argv;

    public:

    class ParseError: public runtime_error
    {
        public:
        ParseError(string const& msg);
    };

    Parser(const int argc, char** const argv);
    Args Parse();
    static void PrintUsage();
};

#endif
