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

#include <alpha/tools/local_poincare4_constraint.hpp>
#include <beta/pcr3bp_reg_basic_objects.hpp>

namespace Pcr3bpProof
{

class SectionPlot4_CE
{
public:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    struct Param
    {
        Pcr3bp::RegBasicObjects<MapT>& basic_objects;
        CapdUtils::LocalCoordinateSystem<MapT> coordsys;
        double span;
        double scale;
        const Leo::Matrix4f& matrix;
        float thickness;
    };

    SectionPlot4_CE(Lyra::Core3d& core_ref, const Param& param)
        : m_core_ref(core_ref)
        , m_param(param)
        , m_linear( param.coordsys.get_origin(), param.coordsys.get_directions_matrix() * param.scale )
        , m_constraint( param.basic_objects.m_hamiltonian_reg2, param.coordsys )
        , m_renderable(
            std::ref(m_composite),
            Leo::RulerSet<2>({
                Leo::Ruler<double>( -param.span, param.span, 21, 1 ),
                Leo::Ruler<double>( -param.span, param.span, 21, 1 ) }),
            std::cref(param.matrix) )
    {
        refresh();
        m_core_ref.register_manifold(&m_renderable);
    }

    ~SectionPlot4_CE() noexcept
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

    CapdUtils::AffineMap<MapT> m_linear;

    LocalPoincare4_Constraint<MapT> m_constraint;

    MapT m_reorder { CapdUtils::ProjectionMap<MapT>::create(4, { 2, 3, 0, 1 }) };

    CapdUtils::CompositeMap<MapT, decltype(m_constraint)&, decltype(m_linear)&, MapT&> m_composite
    {
        std::ref(m_constraint),
        std::ref(m_linear),
        std::ref(m_reorder),
    };

    CapdMapRenderable<decltype(m_composite), 2, 4> m_renderable;
};

}
