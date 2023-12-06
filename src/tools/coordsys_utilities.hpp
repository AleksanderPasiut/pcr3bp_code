///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/local_coordinate_system.hpp>
#include <capd_utils/gauss.hpp>

namespace Carina
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Change coordinate system of the given vector
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class ChangeCoordsys
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static VectorType eval(
        const VectorType& arg,
        const LocalCoordinateSystem<MapT>& src,
        const LocalCoordinateSystem<MapT>& dst)
    {
        const MatrixType A = gaussInverseMatrix<MapT>(src.get_directions_matrix());
        const VectorType b = src.get_origin();
        const MatrixType C = dst.get_directions_matrix();
        const VectorType d = dst.get_origin();

        return (A*C)*arg + A*(d-b);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Generate vector of coordinate systems from another one
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class CoordsysVec
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    template<typename MapU>
    static std::vector<Carina::LocalCoordinateSystem<MapT>> convert( const std::vector<Carina::LocalCoordinateSystem<MapU>>& arg )
    {
        std::vector<Carina::LocalCoordinateSystem<MapT>> ret {};
        ret.reserve(arg.size());

        for (const auto& c : arg)
        {
            ret.emplace_back( Carina::LocalCoordinateSystem<MapT>::convert_from(c));
        }

        return ret;
    }
};


}
