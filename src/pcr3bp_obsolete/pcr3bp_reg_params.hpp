///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"
#include <pcr3bp_basic/setup_parameters.hpp>

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Intermediate parameters for regularized pcr3bp calculations
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Pcr3bpRegParams
{
public:
    using ScalarType = typename MapT::ScalarType;

    ScalarType get_muA() const noexcept
    {
        return m_muA;
    }

    ScalarType get_muB() const noexcept
    {
        return m_muB;
    }

    ScalarType get_xA() const noexcept
    {
        return m_xA;
    }

    ScalarType get_e1() const noexcept
    {
        return m_e1;
    }

    void apply_to(MapT& map) const noexcept
    {
        map.setParameter(0, m_muA);
        map.setParameter(1, m_muB);
        map.setParameter(2, m_xA);
        map.setParameter(3, m_e1);
    }

    static Pcr3bpRegParams<MapT> create(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup)
    {
        switch (mu_index)
        {
            case 1: return Pcr3bpRegParams<MapT>( setup.get_mu(1), setup.get_mu(2), setup.get_x(1), ScalarType(+1.0) );
            case 2: return Pcr3bpRegParams<MapT>( setup.get_mu(2), setup.get_mu(1), setup.get_x(2), ScalarType(-1.0) );
            default: throw std::logic_error("Unsupported i value!");
        }
    }
    
private:
    Pcr3bpRegParams(ScalarType muA, ScalarType muB, ScalarType xA, ScalarType e1)
        : m_muA(muA)
        , m_muB(muB)
        , m_xA(xA)
        , m_e1(e1)
    {}

    ScalarType m_muA;
    ScalarType m_muB;
    ScalarType m_xA;
    ScalarType m_e1;
};

}
