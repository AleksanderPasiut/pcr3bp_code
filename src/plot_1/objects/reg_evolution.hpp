///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>
#include <carina/timemap_wrapper.hpp>
#include <carina/composite_map.hpp>

#include "tools/test_tools.hpp"
#include "tools/direction.hpp"
#include "tools/plotting/solution_curve_interpolation.hpp"

#include "capd_renderable.hpp"

#include "pcr3bp_basic/regularized_system.hpp"
#include "pcr3bp_basic/levi_civita_coordinate_change.hpp"

namespace Ursa
{

struct RegEvolutionParam
{
    Pcr3bp::SetupParameters<RMap> setup;
    Real u0;
    Real v0;
    Real pu0;
    Real pv0;
    Real h;
    Real t;

    size_t point_count;
    size_t point_subcount;

    float line_thickness;
    float point_thickness;

    Leo::Color color;
};

class RegEvolution final
{
public:
    RegEvolution(Lyra::Core2d& core_ref, const RegEvolutionParam& param, Direction direction = Direction::Positive)
        : m_core_ref(core_ref)
        , m_map(create_vector_field(param.setup, direction))
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
    RMap create_vector_field(const Pcr3bp::SetupParameters<RMap>& setup, Direction direction)
    {
        switch (direction)
        {
            case Direction::Positive:
            {
                return Pcr3bp::RegularizedSystem<RMap>::createPositiveVectorField(2, setup, false);
            }
            case Direction::Negative:
            {
                return Pcr3bp::RegularizedSystem<RMap>::createNegativeVectorField(2, setup, false);
            }
            default:
            {
                throw std::logic_error("Unknown direction type!");
            }
        }
    }

    Carina::SolutionCurve<RMap>& get_solution(const RegEvolutionParam& param)
    {
        const RVector U0 = { param.u0, param.v0, param.pu0, param.pv0, param.h };

        m_timemap.set_time(param.t + 0.001);
        m_timemap(U0, m_solution);

        return m_solution;
    }

    Lyra::Core2d& m_core_ref;

    RMap m_map;
    Carina::TimemapWrapper<RMap> m_timemap;
    Carina::SolutionCurve<RMap> m_solution;
    SolutionCurveInterpolation<RMap> m_interpolation;

    CapdSolutionCurvePointRenderable<RMap, 2> m_renderable_points;
    CapdMapRenderable<decltype(m_interpolation), 1, 2> m_renderable_line;
};

class RegEvolutionWithCoordChange final
{
public:
    RegEvolutionWithCoordChange(Lyra::Core2d& core_ref, const RegEvolutionParam& param, Direction direction = Direction::Positive)
        : m_core_ref(core_ref)
        , m_map(create_vector_field(param.setup, direction))
        , m_timemap(m_map, 0.0, 20)
        , m_solution(0.0)
        , m_interpolation(get_solution(param), param.point_count)
        , m_interpolation_with_coord_change( m_interpolation, param.setup )
        , m_renderable_line(
            m_core_ref.get_objects(),
            m_interpolation_with_coord_change,
            Leo::RulerSet<1>({ Leo::Ruler<Real>(0.0, param.t, param.point_count, param.point_subcount) }),
            param.color)
    {
        this->m_renderable_line.fill(param.line_thickness);
    }

private:
    RMap create_vector_field(const Pcr3bp::SetupParameters<RMap>& setup, Direction direction)
    {
        switch (direction)
        {
            case Direction::Positive:
            {
                return Pcr3bp::RegularizedSystem<RMap>::createPositiveVectorField(2, setup, false);
            }
            case Direction::Negative:
            {
                return Pcr3bp::RegularizedSystem<RMap>::createNegativeVectorField(2, setup, false);
            }
            default:
            {
                throw std::logic_error("Unknown direction type!");
            }
        }
    }

    Carina::SolutionCurve<RMap>& get_solution(const RegEvolutionParam& param)
    {
        const RVector U0 = { param.u0, param.v0, param.pu0, param.pv0, param.h };

        m_timemap.set_time(param.t + 0.001);
        m_timemap(U0, m_solution);

        return m_solution;
    }

    Lyra::Core2d& m_core_ref;

    RMap m_map;
    Carina::TimemapWrapper<RMap> m_timemap;
    Carina::SolutionCurve<RMap> m_solution;
    SolutionCurveInterpolation<RMap> m_interpolation;

    class SolutionCurveCoordChange
    {
    public:
        using VectorType = RVector;

        using Node = Carina::Node;

        SolutionCurveCoordChange(
            SolutionCurveInterpolation<RMap>& interpolation_ref,
            const Pcr3bp::SetupParameters<RMap>& setup)
                : m_interpolation_ref(interpolation_ref)
                , m_coord_change( LeviCivitaCoordinateChange<RMap>::create(2, setup, true, true, false) )
        {}

        RVector operator() (RVector arg)
        {
            return m_composite( m_interpolation_ref( arg ) );
        }

    private:
        SolutionCurveInterpolation<RMap>& m_interpolation_ref;
        RMap m_coord_change;

        RMap m_magnify { [](Node, Node in[], int, Node out[], int, Node param[], int) -> void
        {
            out[0] = 10 * in[0];
            out[1] = 10 * in[1];
            
        }, 5, 2, 0};

        Carina::CompositeMap<RMap, RMap&, RMap&> m_composite { std::ref(m_coord_change), std::ref(m_magnify) };

    } m_interpolation_with_coord_change;

    CapdMapRenderable<decltype(m_interpolation_with_coord_change), 1, 2> m_renderable_line;
};

}
