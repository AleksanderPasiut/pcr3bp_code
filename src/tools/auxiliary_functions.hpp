///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/capd/map.hpp>
#include <capd_utils/map_base.hpp>

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief A builder of various functions necessary for the proof
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class AuxiliaryFunctions
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;
    
    static MapT eta(ScalarType L)
    {
        using CapdUtils::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& L = param[0];

            Node& x = in[0];
            Node& y = in[1];

            out[0] = (x - L*y) / (1+L);
            out[1] = (y - L*x) / (1+L);
        };

        MapT map(func, 2, 2, 1);
        map.setParameter(0, L);
		return map;
    }

    static MapT R_Inverse(ScalarType a, ScalarType b)
    {
        using CapdUtils::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& a = param[0];
            Node& b = param[1];

            Node& xi = in[0];
            Node& eta = in[1];

            out[0] = 2*xi / b - 1;
            out[1] = (2*eta - (a+b)) / (b-a);
        };

        MapT map(func, 2, 2, 2);
        map.setParameter(0, a);
        map.setParameter(1, b);
		return map;
    }

    static MapT J()
    {
        using CapdUtils::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            out[0] = in[1];
            out[1] = in[0];
        };

        MapT map(func, 2, 2, 0);
		return map;
    }

    static VectorType S_symmetry(const VectorType& v)
    {
        if (v.dimension() == 4)
        {
            return VectorType{ v[0], -v[1], -v[2], v[3] };
        }
        else
        {
            throw std::logic_error("Unexpected vector dimension!");
        }
    }

    static MapT create_psi0_inverse(const std::array<ScalarType, 2>& d)
    {
        using CapdUtils::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& d1 = param[0];
            Node& d2 = param[1];

            Node& u = in[0];
            Node& pu = in[2];

            out[0] = (u/d1 + pu/d2) / 2;
            out[1] = (u/d1 - pu/d2) / 2;
        };

        MapT map(func, 4, 2, 2);
        map.setParameter(0, d.at(0));
        map.setParameter(1, d.at(1));
		return map;
    }
};

}
