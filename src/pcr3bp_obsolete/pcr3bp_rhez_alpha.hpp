///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/composite_map.hpp>

#include "pcr3bp_rhez_common.hpp"
#include "pcr3bp_rhez_alpha_h.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (u, alpha, h) -> RHEZ
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t mu_index>
class Pcr3bpRhez<MapT, mu_index, CoordType::UAlpha> : public Carina::CompositeMap<MapT, MapT, MapT>
{
public:
    Pcr3bpRhez(const Pcr3bp::SetupParameters<MapT>& setup)
        : Carina::CompositeMap<MapT, MapT, MapT>(
            Pcr3bpRhezCommon<MapT>::create_sc(),
            Pcr3bpRhezCommon<MapT>::create_usch_internal( Pcr3bpRegParams<MapT>::create(mu_index, setup) ))
    {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (v, alpha, h) -> RHEZ
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t mu_index>
class Pcr3bpRhez<MapT, mu_index, CoordType::VAlpha> : public Carina::CompositeMap<MapT, MapT, MapT>
{
public:
    Pcr3bpRhez(const Pcr3bp::SetupParameters<MapT>& setup)
        : Carina::CompositeMap<MapT, MapT, MapT>(
            Pcr3bpRhezCommon<MapT>::create_sc(),
            Pcr3bpRhezCommon<MapT>::create_vsch_internal( Pcr3bpRegParams<MapT>::create(mu_index, setup) ))
    {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! RHEZ -> (u, alpha, h)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t mu_index>
class Pcr3bpRhezInverse<MapT, mu_index, CoordType::UAlpha> : public Carina::CompositeMap<MapT, MapT, MapT, Pcr3bpRhezAlphaH<MapT>>
{
public:
    using ScalarType = typename MapT::ScalarType;

    Pcr3bpRhezInverse(const Pcr3bp::SetupParameters<MapT>& setup)
        : Carina::CompositeMap<MapT, MapT, MapT, Pcr3bpRhezAlphaH<MapT>>(
            Pcr3bpRhezCommon<MapT>::create_extend_u_pu_pv_h( Pcr3bpRegParams<MapT>::create(mu_index, setup) ),
            Pcr3bpRhezCommon<MapT>::create_usch_inverse( Pcr3bpRegParams<MapT>::create(mu_index, setup) ),
            Pcr3bpRhezAlphaH<MapT>() )
    {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! RHEZ -> (v, alpha, h)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t mu_index>
class Pcr3bpRhezInverse<MapT, mu_index, CoordType::VAlpha> : public Carina::CompositeMap<MapT, MapT, MapT, Pcr3bpRhezAlphaH<MapT>>
{
public:
    using ScalarType = typename MapT::ScalarType;

    Pcr3bpRhezInverse(const Pcr3bp::SetupParameters<MapT>& setup)
        : Carina::CompositeMap<MapT, MapT, MapT, Pcr3bpRhezAlphaH<MapT>>(
            Pcr3bpRhezCommon<MapT>::create_extend_v_pu_pv_h( Pcr3bpRegParams<MapT>::create(mu_index, setup) ),
            Pcr3bpRhezCommon<MapT>::create_vsch_inverse( Pcr3bpRegParams<MapT>::create(mu_index, setup) ),
            Pcr3bpRhezAlphaH<MapT>() )
    {}
};

}
