///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/map_base.hpp>

namespace Carina
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Gain map implementation
//! @details Define function G(x) = k * x, where k is some real coefficient.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class GainMap : public MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Constructor
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    GainMap(ScalarType gain, unsigned dimension) : m_gain(gain), m_dimension(dimension)
    {}

    VectorType operator() (const VectorType& vec) override
    {
        this->assert_vector_size(vec, dimension(), "GainMap vec vector size mismatch (1)!");

        return vec * m_gain;
    }

    VectorType operator() (const VectorType& vec, MatrixType& mat) override
    {
        this->assert_vector_size(vec, dimension(), "GainMap vec vector size mismatch (2)!");

        mat = MatrixType::Identity(m_dimension) * m_gain;
        return vec * m_gain;
    }

    unsigned dimension() const noexcept override
    {
        return m_dimension;
    }

    unsigned imageDimension() const noexcept override
    {
        return m_dimension;
    }

    ScalarType get_gain() const
    {
        return m_gain;
    }

private:
    const ScalarType m_gain;
    const unsigned m_dimension;
};

}
