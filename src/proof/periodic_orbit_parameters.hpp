///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/readable_interval.hpp>

#include "tools/types.hpp"
#include "tools/auxiliary_functions.hpp"
#include "pcr3bp_basic/setup_parameters.hpp"

namespace Pcr3bpProof
{

template<typename MapT>
class RegLyapunovCollisionOrbitParameters
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    RegLyapunovCollisionOrbitParameters(const Pcr3bp::SetupParameters<MapT>& setup = Pcr3bp::SetupParameters<MapT>()) : m_setup(setup)
    {
        if (setup.get_mu(1) != ScalarType(1.0) / ScalarType(82.0))
        {
            throw std::logic_error("Unsupported mu value!");
        }
    }

    VectorType get_initial_point() const
    {
        const ScalarType mu = m_setup.get_mu(2);
        return VectorType{ 0.0, 0.0, 0.0, sqrt(mu * 8.0) };
    }

    inline ScalarType get_energy() const;
    inline VectorType get_intermediate_point() const;
    inline VectorType get_image_point() const;

    VectorType get_intermediate_point_neg() const
    {
        return AuxiliaryFunctions<MapT>::S_symmetry( get_intermediate_point() );
    }

private:
    const Pcr3bp::SetupParameters<MapT> m_setup;
};

template<>
inline double RegLyapunovCollisionOrbitParameters<RMap>::get_energy() const
{
    return CapdUtils::ReadableScalar<ScalarType>("bfe6c0fe27bfd02f");
}

template<>
inline Interval RegLyapunovCollisionOrbitParameters<IMap>::get_energy() const
{
    return CapdUtils::ReadableInterval<ScalarType>("bfe6c0fe27bfd196", "bfe6c0fe27bfcec7");
}

template<>
inline RVector RegLyapunovCollisionOrbitParameters<RMap>::get_intermediate_point() const
{
    return VectorType{
        CapdUtils::ReadableScalar<ScalarType>("3ff04f717a03da30"),
        CapdUtils::ReadableScalar<ScalarType>("3fa20bb3ecb359f2"),
        CapdUtils::ReadableScalar<ScalarType>("bff9475c54d58b5f"),
        CapdUtils::ReadableScalar<ScalarType>("bfcc2f1f79908047")
    };
}

template<>
inline IVector RegLyapunovCollisionOrbitParameters<IMap>::get_intermediate_point() const
{
    return VectorType{
        CapdUtils::ReadableInterval<ScalarType>("3ff04f717a03d794", "3ff04f717a03dccb"),
        CapdUtils::ReadableInterval<ScalarType>("3fa20bb3ecb316bb", "3fa20bb3ecb39dd5"),
        CapdUtils::ReadableInterval<ScalarType>("bff9475c54d5922d", "bff9475c54d58475"),
        CapdUtils::ReadableInterval<ScalarType>("bfcc2f1f79909dd1", "bfcc2f1f79906280")
    };
}

template<>
inline RVector RegLyapunovCollisionOrbitParameters<RMap>::get_image_point() const
{
    return VectorType{
        CapdUtils::ReadableScalar<ScalarType>("3fefeb9cf0cc416c"),
        CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
        CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
        CapdUtils::ReadableScalar<ScalarType>("c008cd84a5e72dfd")
    };
}

template<>
inline IVector RegLyapunovCollisionOrbitParameters<IMap>::get_image_point() const
{
    return VectorType{
        CapdUtils::ReadableInterval<ScalarType>("3fefeb9cf0cc3943", "3fefeb9cf0cc4988"),
        CapdUtils::ReadableInterval<ScalarType>("0000000000000000", "0000000000000000"),
        CapdUtils::ReadableInterval<ScalarType>("0000000000000000", "0000000000000000"),
        CapdUtils::ReadableInterval<ScalarType>("c008cd84a5ed1c56", "c008cd84a5e13763")
    };
}

}
