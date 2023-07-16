///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/map_compatibility.hpp>

#include <list>
#include <stdexcept>
#include <sstream>

namespace Carina
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @todo
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, typename MapU, typename... MapV>
class DirectionShifting
{
public:
    static_assert(MapCompatibility<MapT, MapU>::value);

    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static constexpr std::size_t N = sizeof...(MapV) + 1;

    DirectionShifting(ScalarType expansion_factor, MapU map_1, MapV... map_2_args)
        : m_map_1(expansion_factor, map_1)
        , m_map_2(expansion_factor, map_2_args...)
    {}

    std::list<VectorType> eval(const VectorType& first_direction, std::initializer_list<VectorType> points)
    {
        if (points.size() == N)
        {
            return this->eval(first_direction, points.begin(), points.end());
        }
        else
        {
            std::stringstream ss {};
            ss << "Unexpected point list size! ";
            ss << points.size() << ' ' << N;
            throw std::logic_error(ss.str());
        }
    }

    template<typename IteratorType>
    std::list<VectorType> eval(const VectorType& first_direction, IteratorType points_begin, IteratorType points_end)
    {
        if (points_begin != points_end)
        {
            std::list<VectorType> ret = m_map_1.eval(first_direction, points_begin, std::next(points_begin, 1));
            ret.splice( ret.end(), m_map_2.eval(*ret.begin(), ++points_begin, points_end) );
            return ret;
        }
        else
        {
            throw std::logic_error("Empty point list provided!");
        }
    }

private:
    DirectionShifting<MapT, MapU> m_map_1;
    DirectionShifting<MapT, MapV...> m_map_2;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @todo
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, typename MapU>
class DirectionShifting<MapT, MapU>
{
public:
    static_assert(MapCompatibility<MapT, MapU>::value);

    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static constexpr std::size_t N = 1;

    DirectionShifting(ScalarType expansion_factor, MapU map) : m_expansion_factor(expansion_factor), m_map(map)
    {}

    template<typename IteratorType>
    std::list<VectorType> eval(const VectorType& first_direction, IteratorType points_begin, IteratorType points_end)
    {
        if (points_begin != points_end)
        {
            MatrixType der( m_map.imageDimension(), m_map.dimension() );
            m_map( *points_begin, der );
            const VectorType dir = der * first_direction;

            ++points_begin;
            if (points_begin != points_end)
            {
                throw std::logic_error("Unexpected point iterators!");
            }

            return { dir / m_expansion_factor };
        }
        else
        {
            throw std::logic_error("Empty point list provided!");
        }
    }

private:
    ScalarType m_expansion_factor;
    MapU m_map;
};


}
