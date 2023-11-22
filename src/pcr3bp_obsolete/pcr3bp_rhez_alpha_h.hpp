///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/map_base.hpp>
#include <carina/arctan2_map.hpp>

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (u, s, c, h) -> (u, alpha, h)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Pcr3bpRhezAlphaH : Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

public:
    Pcr3bpRhezAlphaH()
    {}

    VectorType operator() (const VectorType& vec, MatrixType& mat)
    {
        this->assert_vector_size(vec, 4, "Pcr3bpRhezAlphaH vec vector size mismatch!");

        MatrixType der(1, 2);
        const VectorType alpha = m_arctan( VectorType{vec[2], vec[1]}, der );

        mat = {
            { 1.0, 0.0, 0.0, 0.0},
            { 0.0, der(1, 2), der(1, 1), 0.0 },
            { 0.0, 0.0, 0.0, 1.0 }
        };

        return VectorType{ vec[0], alpha[0], vec[3] };
    }

    VectorType operator() (const VectorType& vec)
    {
        MatrixType dummy;
        return (*this)(vec, dummy);
    }

    unsigned dimension() const noexcept { return 4; }
    unsigned imageDimension() const noexcept { return 3; }

private:
    Carina::Arctan2<MapT> m_arctan;
};

}
