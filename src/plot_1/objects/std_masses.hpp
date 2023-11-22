///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>
#include "tools/test_tools.hpp"
#include "capd_renderable.hpp"

#include "pcr3bp_obsolete/pcr3bp_setup_values.hpp"

namespace Ursa
{

class StdMasses
{
private:
    Pcr3bpSetupValues<RMap> m_setup;

    CapdVectorRenderable<RVector, 2> m_mass_1;
    CapdVectorRenderable<RVector, 2> m_mass_2;

public:
    StdMasses(
        Lyra::Core2d& core_ref,
        float point_size = 0.01f,
        Leo::Color color_1 = Leo::Color(0.4, 0.5, 0.8),
        Leo::Color color_2 = Leo::Color(0.8, 0.2, 0.2))
            : m_setup()
            , m_mass_1({ m_setup.get_x(1), 0.0 }, color_1)
            , m_mass_2({ m_setup.get_x(2), 0.0 }, color_2)
    {
        m_mass_1.fill(0.01f);
        m_mass_2.fill(0.01f);

        core_ref.register_manifold(&m_mass_1);
        core_ref.register_manifold(&m_mass_2);
    }

    const Pcr3bpSetupValues<RMap>& get_setup() const noexcept
    {
        return m_setup;
    }
};

}
