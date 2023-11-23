///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>
#include "tools/test_tools.hpp"
#include "capd_renderable.hpp"

#include <pcr3bp_basic/setup_parameters.hpp>
#include <proof/homoclinic_orbit_origins_generator.hpp>

namespace Ursa
{

class CoordinateSystemsOrigins
{
public:
    CoordinateSystemsOrigins(
        Lyra::Core2d& core_ref,
        Pcr3bp::SetupParameters<RMap> setup,
        int selected_point= -1,
        float point_size = 0.01f,
        Leo::Color color_1 = Leo::Color(0.1, 0.5, 0.1) )
    {
        HomoclinicOrbitOriginsInitial<RMap> m_homoclinic_orbit_origins_initial {};
        HomoclinicOrbitOriginsGenerator<RMap> generator { m_homoclinic_orbit_origins_initial };

        int index = 0;
        for (auto point : generator.get_points())
        {
            if (selected_point == -1 || index == selected_point)
            {
                m_points.emplace_back(std::ref(core_ref.get_objects()), RVector{ point[0], point[1] }, color_1);
                m_points.rbegin()->fill(point_size);
            }
            ++index;
        }
    }

private:
    std::list<CapdVectorRenderable<RVector, 2>> m_points {};
};

}
