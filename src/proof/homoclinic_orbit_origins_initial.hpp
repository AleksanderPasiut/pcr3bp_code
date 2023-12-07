///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/readable_interval.hpp>

#include "tools/types.hpp"

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Container with homoclinic orbit initial origins
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class HomoclinicOrbitOriginsInitial
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    HomoclinicOrbitOriginsInitial()
    {}

    const std::vector<VectorType>& get_points() const noexcept
    {
        return m_points;
    }

    ScalarType get_total_expansion_factor_pos() const noexcept
    {
        return m_total_expansion_factor_pos;
    }

    ScalarType get_total_expansion_factor_neg() const noexcept
    {
        return m_total_expansion_factor_neg;
    }

private:
    std::vector<VectorType> m_points
    {
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3de90e827829c88e"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3e2a1acc6ddb15f9"),
            CapdUtils::ReadableScalar<ScalarType>("40067d3086e3bf09")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefeb9cf153305e"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3e999e0caa2327ae"),
            CapdUtils::ReadableScalar<ScalarType>("c008cd850bcee9d4")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3e808c01ef984092"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3ec13fce4ef5736c"),
            CapdUtils::ReadableScalar<ScalarType>("40067d30867c8da3")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefeb9e5551793e"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3f30ecc9a62804ac"),
            CapdUtils::ReadableScalar<ScalarType>("c008ce91f4d9e177")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3f15dd90c3f7dc63"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3f56cc4d47f312a1"),
            CapdUtils::ReadableScalar<ScalarType>("40067d2f46afef44")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefef30defcc3f4"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3fc91deae1af399f"),
            CapdUtils::ReadableScalar<ScalarType>("c00bee4375e81c19")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fb304b302711e56"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3ff3f0b1814003a5"),
            CapdUtils::ReadableScalar<ScalarType>("400419ea3b0606b8")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3ff440d7b8e22d64"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3cd1000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3fbec12bfe6c916a")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fb304b302711e56"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bff3f0b1814003a5"),
            CapdUtils::ReadableScalar<ScalarType>("400419ea3b0606b8")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefef30defcc3f4"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bfc91deae1af399f"),
            CapdUtils::ReadableScalar<ScalarType>("c00bee4375e81c19")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3f15dd90c3f7dc63"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bf56cc4d47f312a1"),
            CapdUtils::ReadableScalar<ScalarType>("40067d2f46afef44")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefeb9e5551793e"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bf30ecc9a62804ac"),
            CapdUtils::ReadableScalar<ScalarType>("c008ce91f4d9e177")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3e808c01ef984092"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bec13fce4ef5736c"),
            CapdUtils::ReadableScalar<ScalarType>("40067d30867c8da3")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefeb9cf153305e"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("be999e0caa2327ae"),
            CapdUtils::ReadableScalar<ScalarType>("c008cd850bcee9d4")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3de90e827829c88e"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("be2a1acc6ddb15f9"),
            CapdUtils::ReadableScalar<ScalarType>("40067d3086e3bf09")
        }
    };

    ScalarType m_total_expansion_factor_pos
    {
        CapdUtils::ReadableScalar<ScalarType>("438b3581527bf5b4")
    };

    ScalarType m_total_expansion_factor_neg
    {
        CapdUtils::ReadableScalar<ScalarType>("438b3581527bf5c5")
    };
};



}
