///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <carina/capd/map.hpp>
#include <carina/readable_scalar.hpp>

#include "print_bootstrap.hpp"

namespace Carina
{

template<>
void BootstrapPrint<RMap>::print(std::ostream& ostr, const VectorType& vec)
{
    ostr << "VectorType{\n";

    for (int i = 0; ; ++i)
    {
        ostr << "\tCarina::ReadableScalar<ScalarType>(\"" << std::string(Carina::ReadableScalar<ScalarType>(vec[i])) << "\")";

        if (i+1 < vec.dimension())
        {
            ostr << ",\n";
        }
        else
        {
            ostr << "\n";
            break;
        }
    }
    ostr << "}";
}

template<>
void BootstrapPrint<RMap>::print(std::ostream& ostr, const MatrixType& mat)
{
    const std::pair<unsigned, unsigned> dimension = mat.dimension();
    const std::size_t rows = dimension.first;
    const std::size_t cols = dimension.second;
    const std::size_t count = rows * cols;
    
    ostr << "MatrixType(" << rows << ", " << cols << ", std::array<double, " << count << ">{\n";

    int i = 1;
    for (auto e : mat)
    {
        ostr << "\tCarina::ReadableScalar<ScalarType>(\"" << std::string(Carina::ReadableScalar<double>(e)) << "\")";

        if (i < count)
        {
            ostr << ",\n";
        }
        else
        {
            ostr << "\n";
        }

        ++i;
    }

    ostr << "}.data())";
}

template<>
void BootstrapPrint<IMap>::print(std::ostream& ostr, const VectorType& vec)
{
    throw std::logic_error("Not implemented!");
}

template<>
void BootstrapPrint<IMap>::print(std::ostream& ostr, const MatrixType& mat)
{
    throw std::logic_error("Not implemented!");
}

}
