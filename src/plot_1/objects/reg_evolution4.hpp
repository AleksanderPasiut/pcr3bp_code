///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include <capd_utils/timemap_wrapper.hpp>
#include "tools/test_tools.hpp"
#include "tools/solution_curve_interpolation.hpp"
#include "pcr3bp_basic/regularized_system.hpp"
#include "capd_renderable.hpp"

namespace Pcr3bpProof
{

class RegEvolution4 final
{
public:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    struct Params
    {
        Pcr3bp::SetupParameters<MapT> const & setup;
        VectorType initial_point;
        ScalarType t;
        size_t point_count;
        float thickness;
        bool positive;

        bool operator!= (const Params& arg) const noexcept
        {
            return
                this->initial_point != arg.initial_point ||
                this->t != arg.t ||
                this->point_count != arg.point_count ||
                this->thickness != arg.thickness;
        }
    };

    RegEvolution4(
        Lyra::Core3dObjects& core_objects_ref,
        const Manifold4_Transformation & transformation_ref,
        const Params& params)
            : m_transformation_ref(transformation_ref)
            , m_map(create_map(params.setup, params.positive))
            , m_timemap(m_map, 0.0, 20)
            , m_params(params)
            , m_ruler( Leo::Ruler<ScalarType>(0.0, params.t, params.point_count, 1) )
            , m_renderable(
                core_objects_ref,
                Leo::RulerSet<1>({ m_ruler }))
    {
        refresh();
    }

    void refresh()
    {
        const VectorType U0 = m_params.initial_point;

        m_timemap.set_time(m_params.t);

        CapdUtils::SolutionCurve<MapT> solution(0.0);

        m_timemap(U0, solution);

        Leo::LinearInterpolationNodesList<ScalarType, ScalarType> length_time_nodes;

        ScalarType length = 0;
        for (size_t n = 0; n < m_params.point_count; ++n)
        {
            const ScalarType time = m_ruler.minimum() + n * m_ruler.step();

            const VectorType dU = solution.timeDerivative(time);
            const ScalarType dL = sqrt(1 + capd::vectalg::scalarProduct(dU, dU));

            length_time_nodes.emplace_back( length, time );

            length += dL;
        }

        length_time_nodes.emplace_back( length, m_ruler.maximum() );

        Leo::LinearInterpolationStandard<ScalarType, ScalarType> length_time_interpolation { length_time_nodes };

        auto func = [this, &length_time_interpolation, &length, &solution](const std::array<double, 1>& in) -> std::array<float, 4>
        {
            const ScalarType arg = std::min(in[0] / solution.getRightDomain(), 1.0) * length_time_interpolation.get_max_argument();
            const ScalarType t = length_time_interpolation.evaluate( arg ); 

            const VectorType U = solution(t);

            std::array<double, 4> tmp = { U[2], U[3], U[0], U[1] };
            Lyra::Point4d out = m_transformation_ref.func()(tmp);
            return out;
        };

        m_renderable.fill( func, m_params.thickness );
    }

    const Params& get_params() const noexcept
    {
        return m_params;
    }

private:
    static MapT create_map(Pcr3bp::SetupParameters<MapT> setup, bool positive)
    {
        return positive ?
            Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(2, setup, false) :
            Pcr3bp::RegularizedSystem<MapT>::createNegativeVectorField(2, setup, false);
    }

    const Manifold4_Transformation & m_transformation_ref;

    MapT m_map;
    CapdUtils::TimemapWrapper<MapT> m_timemap;

    Params m_params;

    Leo::Ruler<ScalarType> m_ruler;

    Lyra::Manifold4<1> m_renderable;
};

}
