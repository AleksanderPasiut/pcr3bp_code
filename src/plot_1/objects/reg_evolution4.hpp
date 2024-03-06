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

    struct Param
    {
        Pcr3bp::SetupParameters<MapT> setup;
        VectorType initial_point;
        ScalarType t;
        size_t point_count;
        const Leo::Matrix4f& rotation;
        float thickness;
        bool positive;
    };

    RegEvolution4(Lyra::Core3d& core_ref, const Param& param)
        : m_core_ref(core_ref)
        , m_map(create_map(param.setup, param.positive))
        , m_timemap(m_map, 0.0, 20)
        , m_param(param)
        , m_ruler( Leo::Ruler<ScalarType>(0.0, param.t, param.point_count, 1) )
        , m_renderable(
            Leo::RulerSet<1>({ m_ruler }),
            param.rotation)
    {
        refresh();

        m_core_ref.register_manifold(&m_renderable);
    }

    virtual ~RegEvolution4() noexcept
    {
        m_core_ref.unregister_manifold(&m_renderable);
    }

    void refresh()
    {
        const VectorType U0 = m_param.initial_point;

        m_timemap.set_time(m_param.t);

        CapdUtils::SolutionCurve<MapT> solution(0.0);

        m_timemap(U0, solution);

        Leo::LinearInterpolationNodesList<ScalarType, ScalarType> length_time_nodes;

        ScalarType length = 0;
        for (size_t n = 0; n < m_param.point_count; ++n)
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

            const float u = U[0];
            const float v = U[1];
            const float pu = U[2];
            const float pv = U[3];

            return { pu, pv, u, v };
        };

        m_renderable.fill( func, m_param.thickness );
    }

private:
    static MapT create_map(Pcr3bp::SetupParameters<MapT> setup, bool positive)
    {
        return positive ?
            Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(2, setup, false) :
            Pcr3bp::RegularizedSystem<MapT>::createNegativeVectorField(2, setup, false);
    }

    Lyra::Core3d& m_core_ref;

    MapT m_map;
    CapdUtils::TimemapWrapper<MapT> m_timemap;

    Param m_param;

    Leo::Ruler<ScalarType> m_ruler;

    Lyra::Manifold4<1> m_renderable;
};

}
