///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/local_coordinate_system.hpp>

namespace Ursa
{

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
