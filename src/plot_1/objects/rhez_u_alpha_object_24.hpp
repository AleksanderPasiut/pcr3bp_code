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
//! (u, alpha) -> (pu, pv, u, v)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Rhez_U_Alpha_24 : public Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Rhez_U_Alpha_24(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType h)
        : Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>(
            VectorType{ 0.0, 0.0, h },
            { 0, 1, -1 },
            { 2, 3, 0, 1 },

            std::cref(setup))
    {}
};

class RenderableUAlpha24
{
private:
    Lyra::Core3d& m_core_ref;

    Rhez_U_Alpha_24<RMap> m_map;
    Lyra::Manifold4<2> m_renderable;

    const float m_thickness;

public:
    RenderableUAlpha24(Lyra::Core3d& core_ref, const Rhez_U_Alpha_Param_24& param, const Leo::Matrix4f& rotation)
        : m_core_ref(core_ref)
        , m_map(param.setup, param.h)
        , m_renderable( Leo::RulerSet<2>({ param.u_ruler, param.alpha_ruler }), rotation)
        , m_thickness(param.thickness)
    {
        this->refresh();
        m_core_ref.register_manifold(&m_renderable);
    }

    virtual ~RenderableUAlpha24() noexcept
    {
        m_core_ref.unregister_manifold(&m_renderable);
    }

    void refresh()
    {
        auto func = [this](const std::array<double, 2>& in) -> std::array<float, 4>
        {
            const RVector X = m_map( convert<RVector, 2>(in) );
            return convert<4>(X);
        };

        m_renderable.fill( func, m_thickness );
    }
};

}
