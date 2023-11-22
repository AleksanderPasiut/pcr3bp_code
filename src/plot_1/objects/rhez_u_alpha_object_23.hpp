///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include <carina/enp_map.hpp>
#include "pcr3bp_obsolete/pcr3bp_rhez_alpha.hpp"
#include "capd_renderable.hpp"
#include "rhez_u_alpha_param.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (u, alpha) -> (pu, pv, u)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Rhez_U_Alpha_23 : public Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Rhez_U_Alpha_23(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType h)
        : Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>(
            VectorType{ 0.0, 0.0, h },
            { 0, 1, -1 },
            { 2, 3, 0 },

            std::cref(setup))
    {}
};

class RenderableUAlpha23
{
private:
    Lyra::Core3d& m_core_ref;

    Rhez_U_Alpha_23<RMap> m_map;
    CapdMapRenderable<Rhez_U_Alpha_23<RMap>, 2, 3> m_renderable;

public:
    RenderableUAlpha23(Lyra::Core3d& core_ref, const Rhez_U_Alpha_Param_23& param)
        : m_core_ref(core_ref)
        , m_map(param.setup, param.h)
        , m_renderable(core_ref.get_objects(), m_map, Leo::RulerSet<2>({ param.u_ruler, param.alpha_ruler }), param.color)
    {
        m_renderable.fill(param.thickness);
        // m_core_ref.register_manifold(&m_renderable);
    }

    virtual ~RenderableUAlpha23() noexcept
    {
        // m_core_ref.unregister_manifold(&m_renderable);
    }
};

}
