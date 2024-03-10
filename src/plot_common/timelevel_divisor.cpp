///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "timelevel_divisor.hpp"

namespace Pcr3bpProof
{

TimelevelDivisor::TimelevelDivisor() : m_timepoint()
{}

void TimelevelDivisor::reset()
{
    m_timepoint = std::chrono::system_clock().now();
}

bool TimelevelDivisor::update()
{
    const Timepoint t = std::chrono::system_clock().now();

    if (t - m_timepoint > std::chrono::milliseconds(500))
    {
        m_timepoint = t;
        return true;
    }

    return false;
}

}
