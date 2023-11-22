///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/enp_map.hpp>
#include "pcr3bp_std_poincare.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! py at start -> px at end (with x0 at start fixed)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Pcr3bpStdLyapunovOrbitPy : public Carina::ENP<MapT, Pcr3bpStdPoincare<MapT>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpStdLyapunovOrbitPy(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType x0)
        : Carina::ENP<MapT, Pcr3bpStdPoincare<MapT>>(
            VectorType{ ScalarType(x0), ScalarType(0.0), ScalarType(0.0), ScalarType(0.0) },
            { -1, -1, -1, 0 },
            { 2 },
            std::cref(setup) )
    {}

    ScalarType get_return_time(const VectorType& vec)
    {
        const VectorType v1 = this->extend()(vec);
        return this->internal_map().get_return_time(v1);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (x0, py) at start -> px at end
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Pcr3bpStdLyapunovOrbitX0Py : public Carina::ENP<MapT, Pcr3bpStdPoincare<MapT>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpStdLyapunovOrbitX0Py(const Pcr3bp::SetupParameters<MapT>& setup)
        : Carina::ENP<MapT, Pcr3bpStdPoincare<MapT>>(
            VectorType{ 0.0, 0.0, 0.0, 0.0 },
            { 0, -1, -1, 1 },
            { 2 },
            std::cref(setup) )
    {}
};

}
