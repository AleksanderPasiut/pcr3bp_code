///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <aquila/paramset.hpp>
#include <aquila/paramset_config.hpp>

namespace Pcr3bpProof
{

class CoreInteriorBase
{
public:
    const size_t PARAMSET_CAPACITY = 32;

    CoreInteriorBase() : m_paramset(PARAMSET_CAPACITY)
    {}

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        m_paramset.set(packet_vector);
    }

    double get_param(size_t index)
    {
        return m_paramset.get(index);
    }

private:
    Aquila::Paramset<double> m_paramset;
};

}
