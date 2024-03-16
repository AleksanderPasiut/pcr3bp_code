///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "capd_renderable.hpp"
#include "tools/test_tools.hpp"

#include <capd_utils/local_coordinate_system.hpp>
#include <capd_utils/affine_map.hpp>
#include <capd_utils/enp_map.hpp>

#include <tools/local_poincare4_constraint.hpp>
#include <proof/pcr3bp_reg_basic_objects.hpp>

namespace Pcr3bpProof
{

class SectionPlot4_CE
{
public:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    struct Params
    {
        Pcr3bp::RegBasicObjects<MapT>& basic_objects;
        CapdUtils::LocalCoordinateSystem<MapT> coordsys;
        double span;
        size_t divs;
        const Manifold4_Transformation & transformation_ref;
        float thickness;
    };

    SectionPlot4_CE(Lyra::Core3d& core_ref, const Params& params)
        : m_core_ref(core_ref)
        , m_params(params)
        , m_linear( params.coordsys.get_origin(), params.coordsys.get_directions_matrix() )
        , m_constraint( params.basic_objects.m_hamiltonian_reg2, params.coordsys )
        , m_renderable(
            std::ref(core_ref.get_objects()),
            std::ref(m_composite),
            Leo::RulerSet<2>({
                Leo::Ruler<double>( -params.span, params.span, params.divs, 1 ),
                Leo::Ruler<double>( -params.span, params.span, params.divs, 1 ) }),
            std::cref(params.transformation_ref) )
    {
        refresh();
    }

    void refresh()
    {
        m_renderable.fill(m_params.thickness);
    }

private:
    Lyra::Core3d& m_core_ref;

    Params m_params;

    CapdUtils::AffineMap<MapT> m_linear;

    LocalPoincare4_Constraint<MapT> m_constraint;

    CapdUtils::CompositeMap<MapT, decltype(m_constraint)&, decltype(m_linear)&> m_composite
    {
        std::ref(m_constraint),
        std::ref(m_linear)
    };

    CapdMapRenderable<decltype(m_composite), 2, 4> m_renderable;
};

}
