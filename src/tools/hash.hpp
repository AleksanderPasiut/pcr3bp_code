///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/capd/basic_types.hpp>

namespace std
{

inline size_t combine_hash(size_t x, size_t y)
{
    size_t const tmp = x + 0x9e3779b9 + (y << 6) + (y >> 2);
    return x ^ tmp;
}

template<>
struct hash<CapdUtils::Interval>
{
    size_t operator() (const CapdUtils::Interval& arg) const
    {
        using BoundType = CapdUtils::Interval::BoundType;

        BoundType arr[2];
        arr[0] = arg.leftBound();
        arr[1] = arg.rightBound();

        std::string_view sv(reinterpret_cast<char*>(arr), 2 * sizeof(BoundType));

        return ha(sv);
    }

private:
    std::hash<std::string_view> ha {};
};

template<>
struct hash<CapdUtils::RVector>
{
    size_t operator() (const CapdUtils::RVector& arg) const
    {
        size_t hash = 0;
        for (CapdUtils::Real v : arg)
        {
            hash = combine_hash( hash, ha(v) );
        }
        return hash;
    }

private:
    std::hash<CapdUtils::Real> ha {};
};

template<>
struct hash<CapdUtils::IVector>
{
    size_t operator() (const CapdUtils::IVector& arg) const
    {
        size_t hash = 0;
        for (CapdUtils::Interval v : arg)
        {
            hash = combine_hash( hash, ha(v) );
        }
        return hash;
    }

private:
    std::hash<CapdUtils::Interval> ha {};
};

}
