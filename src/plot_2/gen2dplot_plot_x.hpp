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

#include <carina/local_map.hpp>

namespace Pcr3bpProof
{

class PlotID
{
private:
    RMap m_map;
    Plot22<RMap> m_plot;

    using Params = __decltype(m_plot)::Params;

public:
    PlotID(Lyra::Core2d& core_ref, float thickness)
        : m_map("var: x, y; fun: x, y;")
        , m_plot( core_ref, m_map, get_params(thickness) )
    {}

    static Params get_params(float thickness)
    {
        Params params;
        params.x1_ruler = { -0.5, 0.5, 4, 1 };
        params.x2_ruler = { -0.5, 0.5, 4, 1 };
        params.thickness = thickness;
        params.x1_scale = 1.0;
        params.x2_scale = 1.0;
        return params;
    }
};

template<typename MapT>
class Setup
{
private:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    const Pcr3bp::SetupParameters<MapT> m_setup;

    Pcr3bpRegPoincarePositiveU_XiEta_CE<MapT> m_hlp2;
    Pcr3bpRegPoincareNegativeU_XiEta_CE<MapT> m_hln2;

    using LocalMapPT = CapdUtils::LocalMapObsolete<decltype(m_hlp2)>;
    using LocalMapNT = CapdUtils::LocalMapObsolete<decltype(m_hln2)>;

    LocalMapPT m_local_map_p1;
    LocalMapPT m_local_map_p2;
    LocalMapNT m_local_map_n;

    Pcr3bpRegPoincare2PositiveU_XiEta_CE<MapT> m_flp2;

    using LocalMapFPT = CapdUtils::LocalMapObsolete<decltype(m_flp2)>;

    LocalMapFPT m_local_map_f_p1;
    LocalMapFPT m_local_map_f_p2;

    static VectorType get_fixed_point_1()
    {
        return { 1.0, 0.0 };
    }

    static VectorType get_fixed_point_2()
    {
        return { -1.9979636989004738723, 0.0 };
    }

    static MatrixType get_coord_1()
    {
        return {
            { 0.16843091584405586736, 0.16843091584405586736 },
            { -0.98571346069125609368, +0.98571346069125609368 }
        };
    }

    static MatrixType get_coord_2()
    {
        return {
             { 0.0053537585419491610234, 0.0053537585419491610234 },
             { +0.99998566853204173821, -0.99998566853204173821 }
        };
    }

public:
    Setup()
        : m_setup
        , m_hlp2(m_setup)
        , m_hln2(m_setup)
        , m_local_map_p1( m_hlp2, get_fixed_point_1(), get_fixed_point_2(), get_coord_1(), get_coord_2())
        , m_local_map_p2( m_hlp2, get_fixed_point_2(), get_fixed_point_1(), get_coord_2(), get_coord_1())
        , m_local_map_n( m_hln2, get_fixed_point_1(), get_fixed_point_2(), get_coord_1(), get_coord_2())
        , m_flp2(m_setup)
        , m_local_map_f_p1( m_flp2, get_fixed_point_1(), get_fixed_point_1(), get_coord_1(), get_coord_1())
        , m_local_map_f_p2( m_flp2, get_fixed_point_2(), get_fixed_point_2(), get_coord_2(), get_coord_2())
    {}

    Pcr3bpRegPoincarePositiveU_XiEta_CE<MapT>& get_hlp2()
    {
        return m_hlp2;
    }

    Pcr3bpRegPoincareNegativeU_XiEta_CE<MapT>& get_hln2()
    {
        return m_hln2;
    }

    LocalMapPT& get_local_map_p1()
    {
        return m_local_map_p1;
    }

    LocalMapPT& get_local_map_p2()
    {
        return m_local_map_p2;
    }

    LocalMapNT& get_local_map_n()
    {
        return m_local_map_n;
    }

    Pcr3bpRegPoincare2PositiveU_XiEta_CE<MapT>& get_flp2()
    {
        return m_flp2;
    }

    LocalMapFPT& get_local_map_f_p1()
    {
        return m_local_map_f_p1;
    }

    LocalMapFPT& get_local_map_f_p2()
    {
        return m_local_map_f_p2;
    }
};

class PlotG
{
private:
    Setup<RMap> & m_setup;

    Plot22<__decltype(m_setup.get_local_map_p1())> m_plot;

    using Params = __decltype(m_plot)::Params;

public:
    PlotG(Lyra::Core2d& core_ref, Setup<RMap> & setup, float thickness, float x1_scale, float x2_scale)
        : m_setup(setup)
        , m_plot( core_ref, m_setup.get_local_map_p1(), get_params(thickness, x1_scale, x2_scale) )
    {}

    static Params get_params(float thickness, float x1_scale, float x2_scale)
    {
        Params params;
        params.x1_ruler = { -5.0e-7, 5.0e-7, 10, 4 };
        params.x2_ruler = { -5.0e-7, 5.0e-7, 10, 4 };
        params.thickness = thickness;
        params.x1_scale = x1_scale;
        params.x2_scale = x2_scale;
        return params;
    }
};

class PlotH
{
private:
    Setup<RMap> & m_setup;

    Plot22<__decltype(m_setup.get_local_map_n())> m_plot;

    using Params = __decltype(m_plot)::Params;

public:
    PlotH(Lyra::Core2d& core_ref, Setup<RMap> & setup, float thickness, float x1_scale, float x2_scale)
        : m_setup(setup)
        , m_plot( core_ref, m_setup.get_local_map_n(), get_params(thickness, x1_scale, x2_scale) )
    {}

    static Params get_params(float thickness, float x1_scale, float x2_scale)
    {
        Params params;
        params.x1_ruler = { -5.0e-7, 5.0e-7, 10, 4 };
        params.x2_ruler = { -5.0e-7, 5.0e-7, 10, 4 };
        params.thickness = thickness;
        params.x1_scale = x1_scale;
        params.x2_scale = x2_scale;
        return params;
    }
};

class PlotF
{
private:
    Setup<RMap> & m_setup;

    Plot22<__decltype(m_setup.get_local_map_f_p1())> m_plot;

    using Params = __decltype(m_plot)::Params;

public:
    PlotF(Lyra::Core2d& core_ref, Setup<RMap> & setup, float thickness, float x1_scale, float x2_scale)
        : m_setup(setup)
        , m_plot( core_ref, m_setup.get_local_map_f_p1(), get_params(thickness, x1_scale, x2_scale) )
    {}

    static Params get_params(float thickness, float x1_scale, float x2_scale)
    {
        Params params;
        params.x1_ruler = { -5.0e-9, 5.0e-9, 50, 4 };
        params.x2_ruler = { -1.0e-9, 1.0e-9, 10, 4 };
        params.thickness = thickness;
        params.x1_scale = x1_scale;
        params.x2_scale = x2_scale;
        return params;
    }
};

class X
{
public:
    using MapT = RMap;
    using ScalarType = MapT::ScalarType;
    using VectorType = MapT::VectorType;
    using MatrixType = MapT::MatrixType;

    X(Lyra::Core2d& core_ref, const std::vector<double>& param)
        : m_core_ref(core_ref)
    {
        int idx = 0;
        const size_t plot_type = param[idx++];
        const size_t point_count = param[idx++];
        const double thickness = param[idx++];

        const double span_x1 = pow(10.0, param[idx++] );
        const double span_x2 = pow(10.0, param[idx++] );
        const double scale = pow(10.0, param[idx++] );

        const double shift_x1 = param[idx++] * span_x1;
        const double shift_x2 = param[idx++] * span_x2;

        m_plotid.reset();
        m_plotg.reset();
        m_ploth.reset();
        m_plot_iv.reset();
        m_plot_jv.reset();
        m_plotf.reset();

        switch (plot_type)
        {
            case 0:
            {
                m_plotid = std::make_unique<PlotID>( m_core_ref, thickness );
                break;
            }
            case 1:
            {
                m_plotg = std::make_unique<PlotG>( m_core_ref, std::ref(m_setup), thickness, 2e6, 1e4 );
                break;
            }
            case 2:
            {
                m_ploth = std::make_unique<PlotH>( m_core_ref, std::ref(m_setup), thickness, 1e4, 2e6 );
                break;
            }
            case 3:
            {
                m_plotg = std::make_unique<PlotG>( m_core_ref, std::ref(m_setup), thickness, 1e4, 1e4 );
                m_ploth = std::make_unique<PlotH>( m_core_ref, std::ref(m_setup), thickness, 1e4, 1e4 );
                break;
            }
            case 4:
            {
                const IVector iv = { Interval(-0.5, 0.5), Interval(-0.5, 0.5) };

                PlotIVector::Params params;
                params.color = Leo::Color(0.0, 0.8, 0.4);
                params.thickness = thickness;
                params.x1_scale = 1.0;
                params.x2_scale = 1.0;

                m_plot_iv = std::make_unique<PlotIVector>( m_core_ref, iv, params );
                break;
            }
            case 5:
            {
                const IVector iv = { Interval(-span_x1, span_x1), Interval(-span_x2, span_x2) };
                const IVector jv = m_setup_iv.get_local_map_p1()(iv);

                PlotIVector::Params params;
                params.color = Leo::Color(0.0, 0.8, 0.4);
                params.thickness = thickness;
                params.x1_scale = scale;
                params.x2_scale = scale;
                m_plot_iv = std::make_unique<PlotIVector>( m_core_ref, iv, params );

                params.color = Leo::Color(0.3, 0.3, 0.8);
                m_plot_jv = std::make_unique<PlotIVector>( m_core_ref, jv, params );

                break;
            }
            case 6:
            {
                m_plotf = std::make_unique<PlotF>( m_core_ref, std::ref(m_setup), thickness, 1e11, 1e6 );
                break;
            }
            case 7:
            {
                // scale = 1e10
                // span_x1 = 1e-10.5
                // span_x2 = 1e-13.0

                const IVector iv = { Interval(-span_x1 + shift_x1, span_x1 + shift_x1), Interval(-span_x2 + shift_x2, span_x2 + shift_x2) };
                const IVector jv = m_setup_iv.get_local_map_f_p1()(iv);

                PlotIVector::Params params;
                params.color = Leo::Color(0.0, 0.8, 0.4);
                params.thickness = thickness;
                params.x1_scale = scale;
                params.x2_scale = scale;
                m_plot_iv = std::make_unique<PlotIVector>( m_core_ref, iv, params );

                params.color = Leo::Color(0.3, 0.3, 0.8);
                m_plot_jv = std::make_unique<PlotIVector>( m_core_ref, jv, params );

                break;
            }
            case 8: // full map in point 1
            {
                const IVector iv = { Interval(-span_x1, span_x1), Interval(-span_x2, span_x2) };
                // const IVector jv = m_setup_iv.get_local_map_f_p1()(iv);
                const IVector jvm = m_setup_iv.get_local_map_p1()(iv);
                const IVector jv = m_setup_iv.get_local_map_p2()(jvm);

                PlotIVector::Params params;
                params.color = Leo::Color(0.0, 0.8, 0.4);
                params.thickness = thickness;
                params.x1_scale = scale;
                params.x2_scale = scale;
                m_plot_iv = std::make_unique<PlotIVector>( m_core_ref, iv, params );

                params.color = Leo::Color(0.3, 0.3, 0.8);
                m_plot_jv = std::make_unique<PlotIVector>( m_core_ref, jv, params );

                std::stringstream ss;
                ss.precision(16);
                ss << iv << '\n' << jv << '\n';
                throw std::logic_error(ss.str());
                break;
            }
            case 9: // full map in point 2
            {
                const IVector iv = { Interval(-span_x1, span_x1), Interval(-span_x2, span_x2) };
                // const IVector jv = m_setup_iv.get_local_map_f_p2()(iv); <---------------------------------- why doesn't this work here???
                const IVector jvm = m_setup_iv.get_local_map_p2()(iv);
                const IVector jv = m_setup_iv.get_local_map_p1()(jvm);

                PlotIVector::Params params;
                params.color = Leo::Color(0.0, 0.8, 0.4);
                params.thickness = thickness;
                params.x1_scale = scale;
                params.x2_scale = scale;
                m_plot_iv = std::make_unique<PlotIVector>( m_core_ref, iv, params );

                params.color = Leo::Color(0.3, 0.3, 0.8);
                m_plot_jv = std::make_unique<PlotIVector>( m_core_ref, jv, params );

                std::stringstream ss;
                ss.precision(16);
                ss << iv << '\n' << jv << '\n';
                throw std::logic_error(ss.str());
                break;
            }
            case 10: // half map from point 1 to 2
            {
                const IVector iv = { Interval(-span_x1, span_x1), Interval(-span_x2, span_x2) };
                const IVector jv = m_setup_iv.get_local_map_p1()(iv);

                PlotIVector::Params params;
                params.color = Leo::Color(0.0, 0.8, 0.4);
                params.thickness = thickness;
                params.x1_scale = scale;
                params.x2_scale = scale;
                m_plot_iv = std::make_unique<PlotIVector>( m_core_ref, iv, params );

                params.color = Leo::Color(0.3, 0.3, 0.8);
                m_plot_jv = std::make_unique<PlotIVector>( m_core_ref, jv, params );

                std::stringstream ss;
                ss.precision(16);
                ss << iv << '\n' << jv << '\n';
                throw std::logic_error(ss.str());
                break;
            }
            case 11: // half map from point 2 to 1
            {
                const IVector iv = { Interval(-span_x1, span_x1), Interval(-span_x2, span_x2) };
                const IVector jv = m_setup_iv.get_local_map_p2()(iv);

                PlotIVector::Params params;
                params.color = Leo::Color(0.0, 0.8, 0.4);
                params.thickness = thickness;
                params.x1_scale = scale;
                params.x2_scale = scale;
                m_plot_iv = std::make_unique<PlotIVector>( m_core_ref, iv, params );

                params.color = Leo::Color(0.3, 0.3, 0.8);
                m_plot_jv = std::make_unique<PlotIVector>( m_core_ref, jv, params );

                std::stringstream ss;
                ss.precision(16);
                ss << iv << '\n' << jv << '\n';
                throw std::logic_error(ss.str());
                break;
            }
        }
    }

private:
    Lyra::Core2d& m_core_ref;

    Setup<RMap> m_setup;
    Setup<IMap> m_setup_iv;

    std::unique_ptr<PlotID> m_plotid;
    std::unique_ptr<PlotG> m_plotg;
    std::unique_ptr<PlotH> m_ploth;
    std::unique_ptr<PlotIVector> m_plot_iv;
    std::unique_ptr<PlotIVector> m_plot_jv;
    std::unique_ptr<PlotF> m_plotf;
};

}
