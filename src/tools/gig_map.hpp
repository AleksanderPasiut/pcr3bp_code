///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/map_base.hpp>
#include <carina/map_compatibility.hpp>

namespace Carina
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief GIG map implementation
//! @details Define function g(x) = 1/z f(z*x) for arbitrary function f, where z is input_gain value.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, typename MapU>
class GigMap : public MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static_assert( Carina::MapCompatibility<MapT, MapU>::value );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    GigMap(MapU & map_u, ScalarType input_gain)
        : m_map_u(map_u)
        , m_input_gain(input_gain)
        , m_input_gain_inv( 1.0 / m_input_gain )
    {}

    VectorType operator() (const VectorType& vec) override
    {
        this->assert_vector_size(vec, dimension(), "GigMap vec vector size mismatch (1)!");

        const VectorType arg = vec * m_input_gain;
        const VectorType ret = m_map_u(arg);
        return ret * m_input_gain_inv;
    }

    VectorType operator() (const VectorType& vec, MatrixType& mat) override
    {
        this->assert_vector_size(vec, dimension(), "GigMap vec vector size mismatch (2)!");

        const VectorType arg = vec * m_input_gain;
        const VectorType ret = m_map_u(arg, mat);
        return ret * m_input_gain_inv;
    }

    unsigned dimension() const noexcept override
    {
        return m_map_u.dimension();
    }

    unsigned imageDimension() const noexcept override
    {
        return m_map_u.imageDimension();
    }

    ScalarType get_input_gain() const
    {
        return m_input_gain;
    }

private:
    MapU & m_map_u;

    const ScalarType m_input_gain;
    const ScalarType m_input_gain_inv;
};

}
