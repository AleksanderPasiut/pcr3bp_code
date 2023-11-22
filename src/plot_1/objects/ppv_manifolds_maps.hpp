///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/enp_map.hpp>

#include "pcr3bp_obsolete/pcr3bp_reg_poincare.hpp"
#include "pcr3bp_obsolete/pcr3bp_rhez_alpha.hpp"
#include "pcr3bp_obsolete/pcr3bp_setup_values.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (u, alpha) -> (pu, pv, u)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class RhezUAlpha_Render : public Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    RhezUAlpha_Render(const Pcr3bpSetupValues<MapT>& setup)
        : Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>(
            VectorType{ 0.0, 0.0, setup.get_collision_orbit_h0() },
            { 0, 1, -1 },
            { 2, 3, 0 },

            std::cref(setup))
    {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (u, alpha) -> (u, v, pu, pv, h)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class RhezUAlpha_Render2 : public Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    RhezUAlpha_Render2(const Pcr3bpSetupValues<MapT>& setup)
        : Carina::ENP<MapT, Pcr3bpRhez<MapT, 2, CoordType::UAlpha>>(
            VectorType{ 0.0, 0.0, setup.get_collision_orbit_h0() },
            { 0, 1, -1 },
            { 0, 1, 2, 3, 4 },

            std::cref(setup))
    {}
};

}
