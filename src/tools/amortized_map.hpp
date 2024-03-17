///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "hash.hpp"

#include <capd_utils/capd/map.hpp>

#include <unordered_map>

namespace Pcr3bpProof
{

template<typename MapT>
class AmortizedMap
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    AmortizedMap(MapT& map) : m_map(map)
    {}

    VectorType operator() (const VectorType& arg)
    {
        if (arg.dimension() == m_map.dimension())
        {
            auto it = m_container.find(arg);
            if (it != m_container.end())
            {
                std::cout << "found\n";
                return it->second;
            }

            const VectorType res = m_map(arg);
            m_container[arg] = res;
            return res;
        }
        else
        {
            throw std::logic_error("Unexpected argument dimension!");
        }
    }

private:
    MapT & m_map;

    std::unordered_map<VectorType, VectorType> m_container {};
};

}
