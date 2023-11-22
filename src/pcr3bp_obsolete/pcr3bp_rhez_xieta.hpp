///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/composite_map.hpp>

#include "pcr3bp_rhez_common.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (xi, eta, h) -> RHEZ with v == 0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t mu_index>
class Pcr3bpRhez<MapT, mu_index, CoordType::UXiEta> : public Carina::CompositeMap<MapT, MapT, MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;

    Pcr3bpRhez(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType u0 = ScalarType(-1))
        : Carina::CompositeMap<MapT, MapT, MapT>(
            Pcr3bpRhezCommon<MapT>::create_xe(u0),
            Pcr3bpRhezCommon<MapT>::create_usch_internal( Pcr3bpRegParams<MapT>::create(mu_index, setup) ))
    {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (xi, eta, h) -> RHEZ with u == 0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t mu_index>
class Pcr3bpRhez<MapT, mu_index, CoordType::VXiEta> : public Carina::CompositeMap<MapT, MapT, MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;

    Pcr3bpRhez(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType u0 = ScalarType(-1))
        : Carina::CompositeMap<MapT, MapT, MapT>(
            Pcr3bpRhezCommon<MapT>::create_xe(u0),
            Pcr3bpRhezCommon<MapT>::create_vsch_internal( Pcr3bpRegParams<MapT>::create(mu_index, setup) ))
    {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! RHEZ -> (xi, eta, h) with v == 0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t mu_index>
class Pcr3bpRhezInverse<MapT, mu_index, CoordType::UXiEta> : public Carina::CompositeMap<MapT, MapT, MapT, MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpRhezInverse(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType u0 = ScalarType(-1))
        : Carina::CompositeMap<MapT, MapT, MapT, MapT>(
            Pcr3bpRhezCommon<MapT>::create_extend_u_pu_pv_h( Pcr3bpRegParams<MapT>::create(mu_index, setup) ),
            Pcr3bpRhezCommon<MapT>::create_usch_inverse( Pcr3bpRegParams<MapT>::create(mu_index, setup) ),
            Pcr3bpRhezCommon<MapT>::create_xieta_h(u0))
        , m_u0(u0)
    {}

    VectorType operator() (const VectorType& vec)
    {
        if (vec[0] > m_u0)
        {
            return Carina::CompositeMap<MapT, MapT, MapT, MapT>::operator()(vec);
        }
        else
        {
            return VectorType{ ScalarType(NAN), ScalarType(NAN), ScalarType(NAN) };
        }
    }

    VectorType operator() (const VectorType& vec, MatrixType& der)
    {
        if (vec[0] > m_u0)
        {
            return Carina::CompositeMap<MapT, MapT, MapT, MapT>::operator()(vec, der);
        }
        else
        {
            return VectorType{ ScalarType(NAN), ScalarType(NAN), ScalarType(NAN) };
        }
    }

private:
    const ScalarType m_u0;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! RHEZ -> (xi, eta, h) with u == 0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t mu_index>
class Pcr3bpRhezInverse<MapT, mu_index, CoordType::VXiEta> : public Carina::CompositeMap<MapT, MapT, MapT, MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpRhezInverse(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType u0 = ScalarType(-1))
        : Carina::CompositeMap<MapT, MapT, MapT, MapT>(
            Pcr3bpRhezCommon<MapT>::create_extend_v_pu_pv_h( Pcr3bpRegParams<MapT>::create(mu_index, setup) ),
            Pcr3bpRhezCommon<MapT>::create_vsch_inverse( Pcr3bpRegParams<MapT>::create(mu_index, setup) ),
            Pcr3bpRhezCommon<MapT>::create_xieta_h(u0))
        , m_u0(u0)
    {}

    VectorType operator() (const VectorType& vec)
    {
        if (vec[0] > m_u0)
        {
            return Carina::CompositeMap<MapT, MapT, MapT, MapT>::operator()(vec);
        }
        else
        {
            return VectorType{ ScalarType(NAN), ScalarType(NAN), ScalarType(NAN) };
        }
    }

    VectorType operator() (const VectorType& vec, MatrixType& der)
    {
        if (vec[0] > m_u0)
        {
            return Carina::CompositeMap<MapT, MapT, MapT, MapT>::operator()(vec, der);
        }
        else
        {
            return VectorType{ ScalarType(NAN), ScalarType(NAN), ScalarType(NAN) };
        }
    }

private:
    const ScalarType m_u0;
};

}
