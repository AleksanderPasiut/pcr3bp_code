///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "reg_evolution.hpp"

#include <proof/homoclinic_orbit_origins_generator.hpp>

namespace Ursa
{

class HomoclinicOrbit final
{
public:
    HomoclinicOrbit(Lyra::Core2d& core_ref, double evolution_time, CurveParam curve_param, bool reg = true)
        : m_core_ref(core_ref)
    {
        RegLyapunovCollisionOrbitParameters<RMap> m_parameters { m_setup };

        RegEvolutionParam param;
        param.v0 = 0.0;
        param.h = m_parameters.get_energy();

        param.setup = m_setup;
        param.t = evolution_time;
        param.curve_param = curve_param;
        param.color = Leo::Color(0.0, 0.6, 0.0);

        HomoclinicOrbitOriginsInitial<RMap> m_homoclinic_orbit_origins_initial {};
        HomoclinicOrbitOriginsGenerator<RMap> generator { m_homoclinic_orbit_origins_initial };

        auto points = generator.get_points();

        {
            auto point = *std::next(points.begin(), 8);
            param.u0 = point[0];
            param.pu0 = point[2];
            param.pv0 = point[3];

            if (reg)
            {
                m_segments.emplace_back( std::ref(core_ref), param, Direction::Positive );
            }
            else
            {
                m_segments_std.emplace_back( std::ref(core_ref), param, Direction::Positive );
            }
        }
        {
            auto point = *std::next(points.rbegin(), 8);
            param.u0 = point[0];
            param.pu0 = point[2];
            param.pv0 = point[3];

            if (reg)
            {
                m_segments.emplace_back( std::ref(core_ref), param, Direction::Negative );
            }
            else
            {
                m_segments_std.emplace_back( std::ref(core_ref), param, Direction::Negative );
            }
        }
    }

private:
    Lyra::Core2d& m_core_ref;

    Pcr3bp::SetupParameters<RMap> m_setup {};

    std::list<RegEvolution> m_segments {};
    std::list<RegEvolutionWithCoordChange> m_segments_std {};
};

}
