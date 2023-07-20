///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/capd/map.hpp>
#include <carina/map_base.hpp>

namespace Ursa
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
    
    static MapT Y(ScalarType p)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& p = param[0];

            Node& x = in[0];
            Node& y = in[1];

            out[0] = (x - p*y) / (1+p);
            out[1] = (y - p*x) / (1+p);
        };

        MapT map(func, 2, 2, 1);
        map.setParameter(0, p);
		return map;
    }

    static MapT Y_Inverse(ScalarType p)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& p = param[0];

            Node& x = in[0];
            Node& y = in[1];

            out[0] = (x + p*y) / (1-p);
            out[1] = (y + p*x) / (1-p);
        };

        MapT map(func, 2, 2, 1);
        map.setParameter(0, p);
		return map;
    }

    static MapT R(ScalarType w, ScalarType a, ScalarType b)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& w = param[0];
            Node& a = param[1];
            Node& b = param[2];

            Node& xi = in[0];
            Node& eta = in[1];

            out[0] = w*(xi+1)/2;
            out[1] = ( (a+b)+eta*(b-a) )/2;
        };

        MapT map(func, 2, 2, 3);
        map.setParameter(0, w);
        map.setParameter(1, a);
        map.setParameter(2, b);
		return map;
    }

    static MapT R_T(ScalarType w, ScalarType a, ScalarType b)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& w = param[0];
            Node& a = param[1];
            Node& b = param[2];

            Node& xi = in[1];
            Node& eta = in[0];

            out[1] = w*(xi+1)/2;
            out[0] = ( (a+b)+eta*(b-a) )/2;
        };

        MapT map(func, 2, 2, 3);
        map.setParameter(0, w);
        map.setParameter(1, a);
        map.setParameter(2, b);
		return map;
    }

    static MapT R_Inverse(ScalarType w, ScalarType a, ScalarType b)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& w = param[0];
            Node& a = param[1];
            Node& b = param[2];

            Node& xi = in[0];
            Node& eta = in[1];

            out[0] = 2*xi / w - 1;
            out[1] = (2*eta - (a+b)) / (b-a);
        };

        MapT map(func, 2, 2, 3);
        map.setParameter(0, w);
        map.setParameter(1, a);
        map.setParameter(2, b);
		return map;
    }

    static MapT R_T_Inverse(ScalarType w, ScalarType a, ScalarType b)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& w = param[0];
            Node& a = param[1];
            Node& b = param[2];

            Node& xi = in[1];
            Node& eta = in[0];

            out[1] = 2*xi / w - 1;
            out[0] = (2*eta - (a+b)) / (b-a);
        };

        MapT map(func, 2, 2, 3);
        map.setParameter(0, w);
        map.setParameter(1, a);
        map.setParameter(2, b);
		return map;
    }

    static MapT J()
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            out[0] = in[1];
            out[1] = in[0];
        };

        MapT map(func, 2, 2, 0);
		return map;
    }

    static MapT J2()
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            out[0] = in[0];
            out[1] = -in[1];
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
};

}
