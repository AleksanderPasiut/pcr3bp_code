///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include <carina/enp_map.hpp>
#include "pcr3bp_obsolete/pcr3bp_rhez_xieta.hpp"
#include "capd_renderable.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (xi, eta) -> (pu, pv, u)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Rhez_U_XiEta_23 : public Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UXiEta>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Rhez_U_XiEta_23(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType h)
        : Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UXiEta>>(
            VectorType{ 0.0, 0.0, h },
            { 0, 1, -1 },
            { 2, 3, 0 },

            std::cref(setup), -1.0)
    {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Rhez_U_XiEta_23_Limited
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Rhez_U_XiEta_23_Limited(const Pcr3bp::SetupParameters<MapT>& setup, Real h)
        : m_rhez(setup, h)
    {}

    VectorType operator() (const VectorType& vec)
    {
        const ScalarType xi = vec[0];
        const ScalarType eta = vec[1];

        if ( std::abs(xi*xi + eta*eta - 4) > 0.015)
        {
            return m_rhez(vec);
        }
        else
        {
            return VectorType{ NAN, NAN, NAN };
        }
    }

private:
    Rhez_U_XiEta_23<MapT> m_rhez;
};

struct Rhez_U_Alpha_XiEta_23
{
    Pcr3bp::SetupParameters<RMap> setup;
    Real h;
    Leo::Ruler<Real> xi_ruler;
    Leo::Ruler<Real> eta_ruler;

    float thickness;

    Leo::Color color;
};

class RenderableUXiEta23
{
public:
    RenderableUXiEta23(Lyra::Core3d& core_ref, const Rhez_U_Alpha_XiEta_23& param)
        : m_core_ref(core_ref)
        , m_map(param.setup, param.h)
        , m_renderable(m_map, Leo::RulerSet<2>({ param.xi_ruler, param.eta_ruler }), param.color)
    {
        m_renderable.fill(param.thickness);
        m_core_ref.register_manifold(&m_renderable);
    }

    virtual ~RenderableUXiEta23() noexcept
    {
        m_core_ref.unregister_manifold(&m_renderable);
    }

private:
    Lyra::Core3d& m_core_ref;

    Rhez_U_XiEta_23_Limited<RMap> m_map;
    CapdMapRenderable<Rhez_U_XiEta_23_Limited<RMap>, 2, 3> m_renderable;
};

}
