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
            Carina::ReadableScalar<ScalarType>("3de90e836175052b"),
            Carina::ReadableScalar<ScalarType>("0000000000000000"),
            Carina::ReadableScalar<ScalarType>("3e2a1acd6a283e1f"),
            Carina::ReadableScalar<ScalarType>("40067d3086e3bf09")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3fefeb9cf153305e"),
            Carina::ReadableScalar<ScalarType>("0000000000000000"),
            Carina::ReadableScalar<ScalarType>("3e999e0cbdc16781"),
            Carina::ReadableScalar<ScalarType>("c008cd850bceea24")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3e808c01cf69f720"),
            Carina::ReadableScalar<ScalarType>("0000000000000000"),
            Carina::ReadableScalar<ScalarType>("3ec13fce4f22b9ab"),
            Carina::ReadableScalar<ScalarType>("40067d30867c8e15")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3fefeb9e55517940"),
            Carina::ReadableScalar<ScalarType>("0000000000000000"),
            Carina::ReadableScalar<ScalarType>("3f30ecc9a62dcdb4"),
            Carina::ReadableScalar<ScalarType>("c008ce91f4d9e33c")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3f15dd90c3d60050"),
            Carina::ReadableScalar<ScalarType>("0000000000000000"),
            Carina::ReadableScalar<ScalarType>("3f56cc4d47f331d5"),
            Carina::ReadableScalar<ScalarType>("40067d2f46aff031")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3fefef30defcc3f7"),
            Carina::ReadableScalar<ScalarType>("0000000000000000"),
            Carina::ReadableScalar<ScalarType>("3fc91deae1af46bb"),
            Carina::ReadableScalar<ScalarType>("c00bee4375e81fe9")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3fb304b302710e2e"),
            Carina::ReadableScalar<ScalarType>("0000000000000000"),
            Carina::ReadableScalar<ScalarType>("3ff3f0b1814005ec"),
            Carina::ReadableScalar<ScalarType>("400419ea3b06083c")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3ff440d7b8e22ddf"),
            Carina::ReadableScalar<ScalarType>("0000000000000000"),
            Carina::ReadableScalar<ScalarType>("bcc7000000000000"),
            Carina::ReadableScalar<ScalarType>("3fbec12bfe6c74df")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3fb304b302710e2e"),
            Carina::ReadableScalar<ScalarType>("8000000000000000"),
            Carina::ReadableScalar<ScalarType>("bff3f0b1814005ec"),
            Carina::ReadableScalar<ScalarType>("400419ea3b06083c")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3fefef30defcc3f7"),
            Carina::ReadableScalar<ScalarType>("8000000000000000"),
            Carina::ReadableScalar<ScalarType>("bfc91deae1af46bb"),
            Carina::ReadableScalar<ScalarType>("c00bee4375e81fe9")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3f15dd90c3d60050"),
            Carina::ReadableScalar<ScalarType>("8000000000000000"),
            Carina::ReadableScalar<ScalarType>("bf56cc4d47f331d5"),
            Carina::ReadableScalar<ScalarType>("40067d2f46aff031")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3fefeb9e55517940"),
            Carina::ReadableScalar<ScalarType>("8000000000000000"),
            Carina::ReadableScalar<ScalarType>("bf30ecc9a62dcdb4"),
            Carina::ReadableScalar<ScalarType>("c008ce91f4d9e33c")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3e808c01cf69f720"),
            Carina::ReadableScalar<ScalarType>("8000000000000000"),
            Carina::ReadableScalar<ScalarType>("bec13fce4f22b9ab"),
            Carina::ReadableScalar<ScalarType>("40067d30867c8e15")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3fefeb9cf153305e"),
            Carina::ReadableScalar<ScalarType>("8000000000000000"),
            Carina::ReadableScalar<ScalarType>("be999e0cbdc16781"),
            Carina::ReadableScalar<ScalarType>("c008cd850bceea24")
        },
        VectorType{
            Carina::ReadableScalar<ScalarType>("3de90e836175052b"),
            Carina::ReadableScalar<ScalarType>("8000000000000000"),
            Carina::ReadableScalar<ScalarType>("be2a1acd6a283e1f"),
            Carina::ReadableScalar<ScalarType>("40067d3086e3bf09")
        }
    };

    ScalarType m_total_expansion_factor_pos
    {
        Carina::ReadableScalar<ScalarType>("438b3581527bf58d")
    };

    ScalarType m_total_expansion_factor_neg
    {
        Carina::ReadableScalar<ScalarType>("438b3581527bf590")
    };
};



}
