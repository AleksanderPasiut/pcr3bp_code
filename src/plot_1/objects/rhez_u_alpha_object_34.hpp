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
//! (h, u, alpha) -> (pu, pv, u, h)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Rhez_U_Alpha_34 : public Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Rhez_U_Alpha_34(const Pcr3bp::SetupParameters<MapT>& setup)
        : Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>(
            VectorType{ 0.0, 0.0, 0.0 },
            { 1, 2, 0 },
            { 2, 3, 0, 4 },

            std::cref(setup))
    {}
};

class RenderableUAlpha34
{
private:
    Lyra::Core3d& m_core_ref;

    Rhez_U_Alpha_34<RMap> m_map;
    Lyra::Manifold4<3> m_renderable;

    const float m_thickness;

public:
    RenderableUAlpha34(Lyra::Core3d& core_ref, const Rhez_U_Alpha_Param_34& param, const Leo::Matrix4f& rotation)
        : m_core_ref(core_ref)
        , m_map(param.setup)
        , m_renderable(Leo::RulerSet<3>({ param.h_ruler, param.u_ruler, param.alpha_ruler }), rotation)
        , m_thickness(param.thickness)
    {
        this->refresh();
        m_core_ref.register_manifold(&m_renderable);
    }

    virtual ~RenderableUAlpha34() noexcept
    {
        m_core_ref.unregister_manifold(&m_renderable);
    }

    void refresh()
    {
        auto func = [this](const std::array<double, 3>& in) -> std::array<float, 4>
        {
            const RVector X = m_map( convert<RVector, 3>(in) );
            return convert<4>(X);
        };

        m_renderable.fill( func, m_thickness );
    }
};

}
