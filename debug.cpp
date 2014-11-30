#include "debug.h"

#ifdef WAITFORGDB
#include <unistd.h> // sleep
void WaitForGdb()
{
    const int t = 5; // seconds
    int flag = 1;
    while (flag)
    {
        sleep(t);
    }
}
#endif // WAITFORGDB
