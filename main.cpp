#include <cstdlib>
#include <iostream>

#include "parser.h"
#include "worker.h"

using std::cerr;
using std::endl;

int main(int argc, char** argv)
{
    int exit_code = EXIT_SUCCESS;

    try
    {
        Parser parser(argc, argv);
        Args args = parser.Parse();
        Worker worker(args);
        worker.Run();
    }
    catch (Parser::ParseError& e)
    {
        cerr << e.what() << endl;
        Parser::PrintUsage();
        exit_code = EXIT_FAILURE;
    }

    return exit_code;
}
