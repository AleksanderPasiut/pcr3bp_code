///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>

namespace Carina
{

template<typename MapT>
class BootstrapPrint
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static void print(std::ostream& ostr, const VectorType& vec);
    static void print(std::ostream& ostr, const MatrixType& mat);
};

}
