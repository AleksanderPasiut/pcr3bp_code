///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "setup_parameters.hpp"

namespace Ursa
{
enum class LeviCivitaCoordinateChangeInverseVariant
{
    PositiveU,
    PositiveV,
    NegativeU,
    NegativeV
};


template<typename MapT>
class LeviCivitaCoordinateChange
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Create coordinate change from regularized coordinates to standard coordinates
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create(size_t mu_index, Pcr3bp::SetupParameters<MapT> setup, bool full_change, bool append_h, bool append_t)
    {
        using Carina::Node;

        const ScalarType x0 = setup.get_x(mu_index);

        auto func = [full_change, append_h, append_t](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& x0 = param[0];
            Node& u = in[0];
            Node& v = in[1];

            out[0] = sqr(u) - sqr(v) + x0;
            out[1] = 2*u*v;

            int i = 2;

            if (full_change)
            {
                Node& pu = in[i+0];
                Node& pv = in[i+1];

                Node div = 2 * (sqr(u) + sqr(v));
                out[i+0] = (u*pu-v*pv) / div;
                out[i+1] = (v*pu+u*pv) / div;

                i += 2;
            }

            if (append_h)
            {
                Node& h = in[i];
                out[i] = h;
                i += 1;
            }

            if (append_t)
            {
                Node& t = in[i];
                out[i] = t;
                i += 1;
            }
        };

        const int dimension = 2 + (full_change ? 2 : 0) + (append_h ? 1 : 0) + (append_t ? 1 : 0);
        MapT map(func, dimension, dimension, 1);
        map.setParameter(0, x0);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Create coordinate change from standard to regularized coordinates
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createInverse(
        size_t mu_index,
        Pcr3bp::SetupParameters<MapT> setup,
        bool full_change,
        bool append_h,
        LeviCivitaCoordinateChangeInverseVariant variant = LeviCivitaCoordinateChangeInverseVariant::PositiveU)
    {
        using Carina::Node;

        const ScalarType x0 = setup.get_x(mu_index);

        auto func = [full_change, append_h, variant](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& x0 = param[0];
            Node& x = in[0];
            Node& y = in[1];

            Node R = sqrt( sqr(x-x0) + sqr(y) );

            Node& u = out[0];
            Node& v = out[1];

            switch (variant)
            {
                case LeviCivitaCoordinateChangeInverseVariant::PositiveU:
                {
                    u = sqrt((R + x - x0) / 2 );
                    v = y / (2*u);
                    break;
                }
                case LeviCivitaCoordinateChangeInverseVariant::PositiveV:
                {
                    v = sqrt((R - x + x0) / 2 );
                    u = y / (2*v);
                    break;
                }
                case LeviCivitaCoordinateChangeInverseVariant::NegativeU:
                {
                    u = -sqrt((R + x - x0) / 2 );
                    v = y / (2*u);
                    break;
                }
                case LeviCivitaCoordinateChangeInverseVariant::NegativeV:
                {
                    v = -sqrt((R - x + x0) / 2 );
                    u = y / (2*v);
                    break;
                }
                default:
                {
                    throw std::logic_error("Unknown inverse variant!");
                }
            }

            int i = 2;

            if (full_change)
            {
                Node& px = in[i+0];
                Node& py = in[i+1];

                out[i+0] = 2*(u*px + v*py);
                out[i+1] = 2*(u*py - v*px);

                i += 2;
            }

            if (append_h)
            {
                Node& h = in[i];
                out[i] = h;
                i += 1;
            }
        };

        const int dimension = 2 + (full_change ? 2 : 0) + (append_h ? 1 : 0);
        MapT map(func, dimension, dimension, 1);
        map.setParameter(0, x0);
        return map;
    }

};

}
