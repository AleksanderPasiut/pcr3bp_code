///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "capd_renderable.hpp"
#include "tools/test_tools.hpp"

#include <carina/local_coordinate_system.hpp>
#include <carina/affine_map.hpp>
#include <carina/enp_map.hpp>

namespace Pcr3bpProof
{

class SectionPlot4
{
public:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    struct Param
    {
        Carina::LocalCoordinateSystem<MapT> coordsys;
        double span;
        const Leo::Matrix4f& matrix;
        float thickness;
    };

    SectionPlot4(Lyra::Core3d& core_ref, const Param& param)
        : m_core_ref(core_ref)
        , m_param(param)
        , m_linear( param.coordsys )
        , m_renderable(
            std::ref(m_linear_P0_3),
            Leo::RulerSet<3>({
                Leo::Ruler<double>( -param.span, param.span, 11, 1 ),
                Leo::Ruler<double>( -param.span, param.span, 11, 1 ),
                Leo::Ruler<double>( -param.span, param.span, 11, 1 ) }),
            std::cref(param.matrix) )
    {
        refresh();
        m_core_ref.register_manifold(&m_renderable);
    }

    ~SectionPlot4() noexcept
    {
        m_core_ref.unregister_manifold(&m_renderable);
    }

    void refresh()
    {
        m_renderable.fill(m_param.thickness);
    }

private:
    Lyra::Core3d& m_core_ref;

    Param m_param;

    Carina::AffineMap<MapT> m_linear;

    Carina::ENP<MapT, decltype(m_linear)&> m_linear_P0_3
    {
        VectorType(4),
        { 0, 1, -1, 2 },
        { 2, 3, 0, 1 },
        std::ref(m_linear)
    };

    CapdMapRenderable<decltype(m_linear_P0_3), 3, 4> m_renderable;
};

}
