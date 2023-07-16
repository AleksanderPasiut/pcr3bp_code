///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <stdexcept>

namespace Carina
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Power iteration algorithm implementation
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class PowerIteration
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static VectorType evaluate(const MatrixType& arg, VectorType vec, size_t iterations)
    {
        if (arg.dimension().first != arg.dimension().second)
        {
            throw std::logic_error("Power iteration is implemented for square matrices only!");
        }

        if (arg.dimension().first != vec.dimension())
        {
            throw std::logic_error("Incorrect vector dimension!");
        }

        vec = normalize(vec);

        for (size_t i = 0; i < iterations; ++i)
        {
            vec = arg * vec;
            vec = normalize(vec);
        }

        return vec;
    }

private:
    static VectorType normalize(const VectorType& v)
    {
        return v / v.euclNorm();
    }
};

}
