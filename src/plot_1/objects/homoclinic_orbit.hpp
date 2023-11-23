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
    HomoclinicOrbit(Lyra::Core2d& core_ref, double evolution_time, size_t point_count)
        : m_core_ref(core_ref)
    {
        RegLyapunovCollisionOrbitParameters<RMap> m_parameters { m_setup };


        RegEvolutionParam param;
        param.v0 = 0.0;
        param.h = m_parameters.get_energy();

        param.setup = m_setup;
        param.t = evolution_time;
        param.point_count = point_count;

        param.point_thickness = 0.0f;//5e-3f;
        param.line_thickness = 0.006f;
        param.point_subcount = 10;
        param.color = Leo::Color(1.0, 0.0, 0.0);

        HomoclinicOrbitOriginsInitial<RMap> m_homoclinic_orbit_origins_initial {};
        HomoclinicOrbitOriginsGenerator<RMap> generator { m_homoclinic_orbit_origins_initial };

        auto points = generator.get_points();

        param.u0 = (*points.begin())[0];
        param.pu0 = (*points.begin())[2];
        param.pv0 = (*points.begin())[3];
        m_segments.emplace_back( std::ref(core_ref), param, Direction::Positive);

        param.u0 = (*points.rbegin())[0];
        param.pu0 = (*points.rbegin())[2];
        param.pv0 = (*points.rbegin())[3];
        m_segments.emplace_back( std::ref(core_ref), param, Direction::Negative);
    }

private:
    Lyra::Core2d& m_core_ref;

    Pcr3bp::SetupParameters<RMap> m_setup {};

    std::list<RegEvolution> m_segments {};
};

}
