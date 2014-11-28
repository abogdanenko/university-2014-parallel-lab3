#include "timer.h"
#include "routines.h"

void Timer::Start()
{
    begin = GetWallTime();
}

void Timer::Stop()
{
    end = GetWallTime();
}

double Timer::GetDelta() const
{
    return end - begin;
}


