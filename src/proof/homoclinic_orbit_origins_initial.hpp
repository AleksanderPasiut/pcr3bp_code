///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/readable_interval.hpp>

#include "tools/types.hpp"

namespace Ursa
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
            CapdUtils::ReadableScalar<ScalarType>("3de90e836175052b"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3e2a1acd6a283e1f"),
            CapdUtils::ReadableScalar<ScalarType>("40067d3086e3bf09")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefeb9cf153305e"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3e999e0cbdc16781"),
            CapdUtils::ReadableScalar<ScalarType>("c008cd850bceea24")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3e808c01cf69f720"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3ec13fce4f22b9ab"),
            CapdUtils::ReadableScalar<ScalarType>("40067d30867c8e15")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefeb9e55517940"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3f30ecc9a62dcdb4"),
            CapdUtils::ReadableScalar<ScalarType>("c008ce91f4d9e33c")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3f15dd90c3d60050"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3f56cc4d47f331d5"),
            CapdUtils::ReadableScalar<ScalarType>("40067d2f46aff031")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefef30defcc3f7"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3fc91deae1af46bb"),
            CapdUtils::ReadableScalar<ScalarType>("c00bee4375e81fe9")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fb304b302710e2e"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3ff3f0b1814005ec"),
            CapdUtils::ReadableScalar<ScalarType>("400419ea3b06083c")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3ff440d7b8e22ddf"),
            CapdUtils::ReadableScalar<ScalarType>("0000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bcc7000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("3fbec12bfe6c74df")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fb304b302710e2e"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bff3f0b1814005ec"),
            CapdUtils::ReadableScalar<ScalarType>("400419ea3b06083c")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefef30defcc3f7"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bfc91deae1af46bb"),
            CapdUtils::ReadableScalar<ScalarType>("c00bee4375e81fe9")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3f15dd90c3d60050"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bf56cc4d47f331d5"),
            CapdUtils::ReadableScalar<ScalarType>("40067d2f46aff031")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefeb9e55517940"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bf30ecc9a62dcdb4"),
            CapdUtils::ReadableScalar<ScalarType>("c008ce91f4d9e33c")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3e808c01cf69f720"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("bec13fce4f22b9ab"),
            CapdUtils::ReadableScalar<ScalarType>("40067d30867c8e15")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3fefeb9cf153305e"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("be999e0cbdc16781"),
            CapdUtils::ReadableScalar<ScalarType>("c008cd850bceea24")
        },
        VectorType{
            CapdUtils::ReadableScalar<ScalarType>("3de90e836175052b"),
            CapdUtils::ReadableScalar<ScalarType>("8000000000000000"),
            CapdUtils::ReadableScalar<ScalarType>("be2a1acd6a283e1f"),
            CapdUtils::ReadableScalar<ScalarType>("40067d3086e3bf09")
        }
    };

    ScalarType m_total_expansion_factor_pos
    {
        CapdUtils::ReadableScalar<ScalarType>("438b3581527bf58d")
    };

    ScalarType m_total_expansion_factor_neg
    {
        CapdUtils::ReadableScalar<ScalarType>("438b3581527bf590")
    };
};



}
