///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"
#include <capd_utils/c1_map.hpp>

namespace Pcr3bpProof
{

const Interval I = { -1.0, +1.0 };
const IVector N = { I, I };

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Compute image of set N under specified map and check covering relation conditions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class CoveringRelationCheck
{
public:
    using MapT = IMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    template<typename MapU>
    CoveringRelationCheck(MapU& map)
    {
        assert_with_exception(map.dimension() == 2);
        assert_with_exception(map.imageDimension() == 2);

        CapdUtils::C1_Map<MapT, MapU&> c1_map
        {
            std::ref(map)
        };

        m_img = c1_map(N, m_der);

        const VectorType left = VectorType{ N[0].left(), N[1] };
        m_img_left = c1_map(left);

        const VectorType right = VectorType{ N[0].right(), N[1] };
        m_img_right = c1_map(right);
    }

    bool contraction_condition() const noexcept
    {
        return m_img[1].subset( I );
    }

    bool expansion_condition() const noexcept
    {
        if (m_img_right[0].leftBound() > I.rightBound() && m_img_left[0].rightBound() < I.leftBound())
        {
            return true;
        }

        if (m_img_left[0].leftBound() > I.rightBound() && m_img_right[0].rightBound() < I.leftBound())
        {
            return true;
        }

        return false;
    }

    const VectorType get_img() const noexcept
    {
        return m_img;
    }

    const MatrixType get_der() const noexcept
    {
        return m_der;
    }

    const VectorType get_img_left() const noexcept
    {
        return m_img_left;
    }

    const VectorType get_img_right() const noexcept
    {
        return m_img_right;
    }

private:
    VectorType m_img { VectorType(2) };
    MatrixType m_der { MatrixType(2,2) };
    VectorType m_img_left { VectorType(2) };
    VectorType m_img_right { VectorType(2) };
};

}
