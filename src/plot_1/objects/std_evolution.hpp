///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>
#include <carina/timemap_wrapper.hpp>
#include "tools/test_tools.hpp"
#include "pcr3bp_basic/standard_system.hpp"
#include "capd_renderable.hpp"

#include "tools/plotting/solution_curve_interpolation.hpp"

namespace Ursa
{

struct StdEvolutionParam
{
    Pcr3bp::SetupParameters<RMap> setup;
    Real x0;
    Real y0;
    Real px0;
    Real py0;
    Real t;

    size_t point_count;
    size_t point_subcount;

    float line_thickness;
    float point_thickness;

    Leo::Color color;
};

class StdEvolution final
{
private:
    Lyra::Core2d& m_core_ref;

    RMap m_map;
    Carina::TimemapWrapper<RMap> m_timemap;
    Carina::SolutionCurve<RMap> m_solution;
    SolutionCurveInterpolation<RMap> m_interpolation;

    CapdSolutionCurvePointRenderable<RMap, 2> m_renderable_points;
    CapdMapRenderable<decltype(m_interpolation), 1, 2> m_renderable_line;

public:
    StdEvolution(Lyra::Core2d& core_ref, const StdEvolutionParam& param)
        : m_core_ref(core_ref)
        , m_map(Pcr3bp::StandardSystem<RMap>::createPositiveVectorField(param.setup, false))
        , m_timemap(m_map, 0.0, 20)
        , m_solution(0.0)
        , m_interpolation(get_solution(param), param.point_count)
        , m_renderable_points(
            core_ref.get_objects(),
            m_solution,
            0.0,
            param.t,
            param.point_count,
            param.color)
        , m_renderable_line(
            core_ref.get_objects(),
            m_interpolation,
            Leo::RulerSet<1>({ Leo::Ruler<Real>(0.0, param.t, param.point_count, param.point_subcount) }),
            param.color)
    {
        this->m_renderable_points.fill(param.point_thickness);
        this->m_renderable_line.fill(param.line_thickness);
    }

private:
    Carina::SolutionCurve<RMap>& get_solution(const StdEvolutionParam& param)
    {
        const RVector X0 = { param.x0, param.y0, param.px0, param.py0 };

        m_timemap.set_time(param.t + 0.001);
        m_timemap(X0, m_solution);

        return m_solution;
    }
};

}
