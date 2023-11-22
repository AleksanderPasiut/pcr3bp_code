///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <leo/ruler.hpp>
#include <leo/misc/color.hpp>
#include "tools/test_tools.hpp"
#include "pcr3bp_obsolete/pcr3bp_setup_values.hpp"

namespace Ursa
{

struct Rhez_U_Alpha_Param_23
{
    Pcr3bpSetupValues<RMap> setup;
    Real h;
    Leo::Ruler<Real> u_ruler;
    Leo::Ruler<Real> alpha_ruler;

    float thickness;

    Leo::Color color;
};

struct Rhez_U_Alpha_Param_24
{
    Pcr3bpSetupValues<RMap> setup;
    Real h;
    Leo::Ruler<Real> u_ruler;
    Leo::Ruler<Real> alpha_ruler;

    float thickness;
};

struct Rhez_U_Alpha_Param_34
{
    Pcr3bpSetupValues<RMap> setup;
    Leo::Ruler<Real> h_ruler;
    Leo::Ruler<Real> u_ruler;
    Leo::Ruler<Real> alpha_ruler;

    float thickness;
};

}
