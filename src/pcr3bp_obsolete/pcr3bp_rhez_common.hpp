///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <pcr3bp_basic/setup_parameters.hpp>
#include "pcr3bp_reg_params.hpp"

namespace Ursa
{

enum class CoordType
{
    UAlpha,
    VAlpha,
    UXiEta,
    VXiEta
};

template<typename MapT, size_t mu_index, CoordType coord_type>
class Pcr3bpRhez
{};

template<typename MapT, size_t mu_index, CoordType coord_type>
class Pcr3bpRhezInverse
{};

class Pcr3bpRhezBase
{
protected:
    using Node = Carina::Node;

    static Node get_ru(Node& u, Node& h, Node& muA, Node& muB, Node& e1, Node& xA)
    {
        return 2*sqrt( 2*sqr(u)*h + 2*muA + 2*muB*sqr(u) / sqrt(sqr(sqr(u) + e1)) + sqr(u)*sqr(sqr(u) + xA));
    }

    static Node get_rv(Node& v, Node& h, Node& muA, Node& muB, Node& e1, Node& xA)
    {
        return 2*sqrt( 2*sqr(v)*h + 2*muA + 2*muB*sqr(v) / sqrt(sqr(sqr(v) - e1)) + sqr(v)*sqr(sqr(v) - xA));
    }

    static Node get_pus(Node& u, Node& xA)
    {
        return 2*u*(sqr(u) + xA);
    }

    static Node get_pvs(Node& v, Node& xA)
    {
        return 2*v*(sqr(v) - xA);
    }
};

template<typename MapT>
class Pcr3bpRhezCommon : public Pcr3bpRhezBase
{
public:
    using ScalarType = typename MapT::ScalarType;

    using Base = Pcr3bpRhezBase;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! (u, alpha, h) -> (u, s, c, h)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create_sc()
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& u = in[0];
            Node& alpha = in[1];
            Node& h = in[2];

            out[0] = u;
            out[1] = sin(alpha);
            out[2] = cos(alpha);
            out[3] = h;
        };

        MapT map(func, 3, 4, 0);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! (xi, eta, h) -> (u, s, c, h)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create_xe(ScalarType u0)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& u0 = param[0];
            Node& xi = in[0];
            Node& eta = in[1];
            Node& h = in[2];

            Node u = sqrt( sqr(xi) + sqr(eta) );

            out[0] = u + u0;
            out[1] = eta / u;
            out[2] = xi / u;
            out[3] = h;
        };

        MapT map(func, 3, 4, 1);
        map.setParameter(0, u0);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! (u, s, c, h) -> (u, 0, pu, pv, h)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create_usch_internal(const Pcr3bpRegParams<MapT>& params)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& muA = param[0];
            Node& muB = param[1];
            Node& xA = param[2];
            Node& e1 = param[3];

            Node& u = in[0];
            Node& s = in[1];
            Node& c = in[2];
            Node& h = in[3];

            Node ru  = get_ru(u, h, muA, muB, e1, xA);
            Node pus = get_pus(u, xA);

            out[0] = u;
            out[1] = Node(0);
            out[2] = ru * s;
            out[3] = ru * c + pus;
            out[4] = h;
        };

        MapT map(func, 4, 5, 4);
        params.apply_to(map);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! (v, s, c, h) -> (0, v, pu, pv, h)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create_vsch_internal(const Pcr3bpRegParams<MapT>& params)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& muA = param[0];
            Node& muB = param[1];
            Node& xA = param[2];
            Node& e1 = param[3];

            Node& v = in[0];
            Node& s = in[1];
            Node& c = in[2];
            Node& h = in[3];

            Node rv  = get_rv(v, h, muA, muB, e1, xA);
            Node pvs = get_pvs(v, xA);

            out[0] = Node(0);
            out[1] = v;
            out[2] = rv * c - pvs;
            out[3] = rv * s;
            out[4] = h;
        };

        MapT map(func, 4, 5, 4);
        params.apply_to(map);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! (u, 0, pu, pv, h) -> (u, pu, pv, ru, pus, h)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create_extend_u_pu_pv_h(const Pcr3bpRegParams<MapT>& params)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& muA = param[0];
            Node& muB = param[1];
            Node& xA = param[2];
            Node& e1 = param[3];

            Node& u = in[0];
            Node& v = in[1];
            Node& pu = in[2];
            Node& pv = in[3];
            Node& h = in[4];

            Node ru  = get_ru(u, h, muA, muB, e1, xA);
            Node pus = get_pus(u, xA);

            out[0] = u;
            out[1] = pu;
            out[2] = pv;
            out[3] = ru;
            out[4] = pus;
            out[5] = h;
        };

        MapT map(func, 5, 6, 4);
        params.apply_to(map);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! (0, v, pu, pv, h) -> (v, pu, pv, rv, pvs, h)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create_extend_v_pu_pv_h(const Pcr3bpRegParams<MapT>& params)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& muA = param[0];
            Node& muB = param[1];
            Node& xA = param[2];
            Node& e1 = param[3];

            Node& u = in[0];
            Node& v = in[1];
            Node& pu = in[2];
            Node& pv = in[3];
            Node& h = in[4];

            Node rv  = get_rv(v, h, muA, muB, e1, xA);
            Node pvs = get_pvs(v, xA);

            out[0] = v;
            out[1] = pu;
            out[2] = pv;
            out[3] = rv;
            out[4] = pvs;
            out[5] = h;
        };

        MapT map(func, 5, 6, 4);
        params.apply_to(map);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! (u, pu, pv, ru, pus, h) -> (u, s, c, h)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create_usch_inverse(const Pcr3bpRegParams<MapT>& params)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& u = in[0];
            Node& pu = in[1];
            Node& pv = in[2];
            Node& ru = in[3];
            Node& pus = in[4];
            Node& h = in[5];

            out[0] = u;
            out[1] = pu / ru;
            out[2] = (pv - pus) / ru;
            out[3] = h;
        };

        MapT map(func, 6, 4, 0);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! (v, pu, pv, rv, pvs, h) -> (v, s, c, h)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create_vsch_inverse(const Pcr3bpRegParams<MapT>& params)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& v = in[0];
            Node& pu = in[1];
            Node& pv = in[2];
            Node& rv = in[3];
            Node& pvs = in[4];
            Node& h = in[5];

            out[0] = v;
            out[1] = pv / rv;
            out[2] = (pu + pvs) / rv;
            out[3] = h;
        };

        MapT map(func, 6, 4, 0);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! (u, s, c, h) -> (xi, eta, h)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT create_xieta_h(ScalarType u0)
    {
        using Carina::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& u0 = param[0];
            Node& u = in[0];
            Node& s = in[1];
            Node& c = in[2];
            Node& h = in[3];

            out[0] = c * (u - u0);
            out[1] = s * (u - u0);
            out[2] = h;
        };

        MapT map(func, 4, 3, 1);
        map.setParameter(0, u0);
        return map;
    }
};

}
