///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>
#include "tools/test_tools.hpp"
#include "capd_renderable.hpp"

#include <pcr3bp_basic/setup_parameters.hpp>
#include <proof/homoclinic_orbit_origins_generator.hpp>
#include <pcr3bp_basic/levi_civita_coordinate_change.hpp>

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
        bool reg = true,
        Leo::Color color_1 = Leo::Color(0.0, 0.0, 0.0))
    {
        HomoclinicOrbitOriginsInitial<RMap> m_homoclinic_orbit_origins_initial {};
        HomoclinicOrbitOriginsGenerator<RMap> generator { m_homoclinic_orbit_origins_initial };

        if (reg)
        {
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
        else
        {
            using Carina::Node;
            RMap coord_change( LeviCivitaCoordinateChange<RMap>::create(2, setup, true, false, false) );
            RMap magnify( [](Node, Node in[], int, Node out[], int, Node param[], int) -> void
            {
                out[0] = 10 * in[0];
                out[1] = 10 * in[1];
                
            }, 4, 2, 0);

            Carina::CompositeMap<RMap, RMap&, RMap&> composite { std::ref(coord_change), std::ref(magnify) };

            for (auto point : generator.get_points())
            {
                RVector point_std = composite(point);
                m_points.emplace_back(std::ref(core_ref.get_objects()), point_std, color_1);
                m_points.rbegin()->fill(point_size);
            }
        }
    }

private:
    std::list<CapdVectorRenderable<RVector, 2>> m_points {};
};

}
