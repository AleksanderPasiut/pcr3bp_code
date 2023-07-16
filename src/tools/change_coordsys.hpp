///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/local_coordinate_system.hpp>
#include <carina/gauss.hpp>

namespace Carina
{

template<typename MapT>
class ChangeCoordsys
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static VectorType eval(
        const VectorType& arg,
        const Carina::LocalCoordinateSystem<MapT>& src,
        const Carina::LocalCoordinateSystem<MapT>& dst)
    {
        const MatrixType A = gaussInverseMatrix<MapT>(src.get_directions_matrix());
        const VectorType b = src.get_origin();
        const MatrixType C = dst.get_directions_matrix();
        const VectorType d = dst.get_origin();

        return (A*C)*arg + A*(d-b);
    }
};

}
