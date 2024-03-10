///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

namespace Pcr3bpProof
{

class EdgeDetector
{
public:
    void reset(bool arg)
    {
        m_previous_value = arg;
    }

    bool update(bool arg)
    {
        bool const ret = (arg != m_previous_value);
        m_previous_value = arg;
        return ret;
    }

private:
    bool m_previous_value {};
};


}
