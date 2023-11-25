///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>
#include "tools/test_tools.hpp"
#include "capd_renderable.hpp"

#include <pcr3bp_basic/setup_parameters.hpp>

namespace Ursa
{

class StdMasses
{
public:
    StdMasses(
        Lyra::Core2d& core_ref,
        Pcr3bp::SetupParameters<RMap> setup,
        float point_size = 0.01f,
        Leo::Color color_1 = Leo::Color(0.4, 0.5, 0.8),
        Leo::Color color_2 = Leo::Color(0.8, 0.2, 0.2))
            : m_mass_1(core_ref.get_objects(), { setup.get_x(1), 0.0 }, color_1)
            , m_mass_2(core_ref.get_objects(), { setup.get_x(2), 0.0 }, color_2)
    {
        m_mass_1.fill(point_size);
        m_mass_2.fill(point_size);
    }

private:
    CapdVectorRenderable<RVector, 2> m_mass_1;
    CapdVectorRenderable<RVector, 2> m_mass_2;
};

}
