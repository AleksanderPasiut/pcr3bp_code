///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

namespace Ursa
{

template<typename MapT>
class ScaleMap
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    ScaleMap(const VectorType& scale) : m_scale(scale)
    {}

    VectorType operator() (const VectorType& arg)
    {
        return { arg[0] * m_scale[0], arg[1] * m_scale[1] };
    }

private:
    VectorType m_scale;
};

}
