///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>

#include "plot22.hpp"
#include "plot_ivector.hpp"

#include "pcr3bp_obsolete/pcr3bp_reg_poincare_negative_u_xieta.hpp"
#include "pcr3bp_obsolete/pcr3bp_reg_poincare_positive_u_xieta.hpp"

#include "pcr3bp_basic/pcr3bp_reg_poincare2_negative_u_xieta.hpp"
#include "pcr3bp_basic/pcr3bp_reg_poincare2_positive_u_xieta.hpp"

#include <carina/identity_map.hpp>
#include <carina/local_map.hpp>

#include "gen2dplot_plot.scale.hpp"

namespace Pcr3bpProof
{

// template<template<typename> typename Func>
// class FI
// {
// public:
//     using MapT = RMap;
//     using ScalarType = typename MapT::ScalarType;
//     using VectorType = typename MapT::VectorType;
//     using MatrixType = typename MapT::MatrixType;
//
//     FI( Func<MapT>& func ) : m_func(func)
//     {}
//
//
//
// private:
//     Func<MapT>& m_func;
//
// };

template<typename MapT, template<typename> typename GMapT>
class W
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    W(
        GMapT<MapT>& gen_map,
        const VectorType& origin,
        const VectorType& origin_image,
        const VectorType& pos_dir,
        const VectorType& neg_dir)
            : m_gen_map(gen_map)
            , m_local(
                m_gen_map,
                origin,
                origin_image,
                c2x2mr(pos_dir, neg_dir),
                c2x2mr(pos_dir, neg_dir) )
    {}

    VectorType operator() (const VectorType& arg)
    {
        return m_local(arg);
    }

private:
    static MatrixType c2x2mr(const VectorType& v1, const VectorType& v2)
    {
        if (v1.dimension() == 2)
        {
            if (v2.dimension() == 2)
            {
                return MatrixType{ { v1[0], v2[0] }, { v1[1], v2[1] } };
            }
            else
            {
                throw std::invalid_argument("Vector v1 dimension mismatch!");
            }
        }
        else
        {
            throw std::invalid_argument("Vector v1 dimension mismatch!");
        }
    }


    GMapT<MapT>& m_gen_map;
    Carina::LocalMapObsolete<GMapT<MapT>> m_local;
};

template<typename MapT>
class V
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    V()
        : m_setup()
        , m_xieta(m_setup)
        , m_w(
            m_xieta,
            VectorType{ 1, 0 },
            VectorType{ 1, 0 },
            get_ppv_positive_unstable_dir<VectorType>(),
            get_ppv_negative_unstable_dir<VectorType>())
    {}

    VectorType operator() (const VectorType& arg)
    {
        return m_w(arg);
    }

private:
    Pcr3bp::SetupParameters<MapT> m_setup;
    Pcr3bpRegPoincare2PositiveU_XiEta_CE<MapT> m_xieta;
    W<MapT, Pcr3bpRegPoincare2PositiveU_XiEta_CE> m_w;
};

template<typename MapT>
class U
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    U(const VectorType& scale)
        : m_v()
        , m_scale(scale)
    {}

    VectorType operator() (const VectorType& arg)
    {
        return m_scale(m_v(arg));
    }

private:
    V<MapT> m_v;
    ScaleMap<MapT> m_scale;
};

class E
{
public:
    E(Real x1, Real x2)
        : m_vi()
        , m_cborder{ Interval{-x1, +x1}, Interval{-x2, +x2} }
        , m_lborder{ Interval{-x1, -x1}, Interval{-x2, +x2} }
        , m_rborder{ Interval{+x1, +x1}, Interval{-x2, +x2} }
        , m_cimg( m_vi(m_cborder) )
        , m_limg( m_vi(m_lborder) )
        , m_rimg( m_vi(m_rborder) )
    {
        std::cout << "cimg\n" << m_cimg << '\n';
        std::cout << "limg\n" << m_limg << '\n';
        std::cout << "rimg\n" << m_rimg << '\n';
    }

    const IVector& get_cborder() const noexcept
    {
        return m_cborder;
    }

    const IVector& get_lborder() const noexcept
    {
        return m_lborder;
    }

    const IVector& get_rborder() const noexcept
    {
        return m_rborder;
    }

    const IVector& get_cimg() const noexcept
    {
        return m_cimg;
    }

    const IVector& get_limg() const noexcept
    {
        return m_limg;
    }

    const IVector& get_rimg() const noexcept
    {
        return m_rimg;
    }

private:
    V<IMap> m_vi;

    IVector m_cborder;
    IVector m_lborder;
    IVector m_rborder;
    IVector m_cimg;
    IVector m_limg;
    IVector m_rimg;
};

class G
{
public:
    G(Lyra::Core2d& core_ref, IVector cborder, Real scale, float thickness)
        : m_core_ref(core_ref)
        , m_u( RVector{ scale, scale } )
        , m_renderable_u(m_u, cborder, 5, 10, Leo::Color(1.0, 0.0, 0.0))
    {
        m_renderable_u.fill( thickness );
        m_core_ref.register_manifold(&m_renderable_u);
    }

    virtual ~G() noexcept
    {
        m_core_ref.unregister_manifold(&m_renderable_u);
    }

private:
    Lyra::Core2d& m_core_ref;

    U<RMap> m_u;
    CapdMapIntervalRenderable<decltype(m_u), 2, 2> m_renderable_u;
};

class X
{
public:
    X(Lyra::Core2d& core_ref, const std::vector<double>& param)
        : m_core_ref(core_ref)
        , m_scale( pow(10, -param[5]) )
        , m_x1( pow(10, param[3]) )
        , m_x2( pow(10, param[4]) )
        , m_e(m_x1, m_x2)
        , m_ivs()
        , m_g( core_ref, m_e.get_cborder(), m_scale, param[2] )
    {
        PlotIVector::Params params;
        params.thickness = param[2];
        params.x1_scale = m_scale;
        params.x2_scale = m_scale;

        params.color = Leo::Color(1.0, 0.0, 1.0); // pink
        m_ivs.emplace_back(std::ref(core_ref), std::cref(m_e.get_cborder()), params);

        params.color = Leo::Color(0.0, 1.0, 1.0); // light blue
        m_ivs.emplace_back(std::ref(core_ref), std::cref(m_e.get_cimg()), params);

        params.color = Leo::Color(0.0, 1.0, 0.0); // green
        m_ivs.emplace_back(std::ref(core_ref), std::cref(m_e.get_limg()), params);

        params.color = Leo::Color(1.0, 1.0, 0.0); // yellow
        m_ivs.emplace_back(std::ref(core_ref), std::cref(m_e.get_rimg()), params);
    }

private:
    Lyra::Core2d& m_core_ref;

    const Real m_scale;
    const Real m_x1;
    const Real m_x2;

    E m_e;

    std::list<PlotIVector> m_ivs;

    G m_g;
};

}
