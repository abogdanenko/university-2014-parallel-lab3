#include <iostream>
#include <sstream>
#include <cctype>
#include <unistd.h>

#include "parser.h"
#include "routines.h"

using std::cout;
using std::endl;
using std::ostringstream;
using std::hex;

Parser::ParseError::ParseError(const string& msg):
    runtime_error(msg)
{

}

Parser::Parser(const int argc, char** const argv):
    argc(argc),
    argv(argv)
{

}

void Parser::PrintUsage()
{
    cout << "Usage: lab3 "
            "[-h] "
            "[-v] "
            "[-n size] "
            "[-e epsilon_max] "
            "[-i iteration_count_max] "
            "[-o omega] "
            "[-a array_output_file] "
            "[-s stats_file]"
        << endl;
}

Args Parser::Parse()
{
    Args result;
    ostringstream oss;
    int c; // option character
    while ((c = getopt(argc, argv, ":hvn:e:o:i:a:s:")) != -1)
    {
        switch(c)
        {
            case 'h':
                result.usage_flag = true;
                break;
            case 'v':
                result.verbose_flag = true;
                break;
            case 'n':
                result.n = string_to_number<unsigned>(optarg);
                break;
            case 'e':
                result.epsilon_max = string_to_number<double>(optarg);
                break;
            case 'o':
                result.omega = string_to_number<double>(optarg);
                break;
            case 'i':
                result.iteration_count_max = string_to_number<unsigned>(optarg);
                break;
            case 'a':
                result.array_filename = optarg;
                break;
            case 's':
                result.stats_filename = optarg;
                break;
            case ':':
                oss << "Option -" << char(optopt) << " requires an argument.";
                throw ParseError(oss.str());
            case '?':
                if (isprint(optopt))
                {
                    oss << "Unknown option `-" << char(optopt) << "'.";
                }
                else
                {
                    oss << "Unknown option character `\\x" << hex << optopt <<
                        "'.";
                }
                throw ParseError(oss.str());
            default:
                throw ParseError("An error occured while parsing options.");
        }
    }

    if (optind < argc)
    {
        throw ParseError("Extra non-option arguments found");
    }

    return result;
}
