///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_3d.hpp>
#include <lyra/core3d.hpp>

#include "plot_common/window_properties.hpp"
#include "plot_common/core_interior_base.hpp"

#include "objects/rhez_u_alpha_object_23.hpp"
#include "objects/rhez_u_xieta_object_23.hpp"
#include "objects/ppv_manifolds.hpp"

#include "objects/interp.hpp"
#include "objects/unstable_ppv_manifold.hpp"

#include "objects/intersection_renderable.hpp"

#include "rhez_u_plot.v3.hpp"
#include "rhez_u_plot.manifold_intersections.hpp"

namespace Ursa
{

class LineManifolds
{
public:
    LineManifolds(Lyra::Core3d& core_ref, const Pcr3bpSetupValues<RMap>& setup)
        : m_core_ref(core_ref)
        , m_setup(setup)
        , m_unstable_pos()
        , m_unstable_neg()
        , m_intersections()
    {}

    void refresh()
    {
        m_unstable_pos = std::make_unique<UnstableManifoldPos>(
            std::ref(m_core_ref),
            std::cref(m_setup),
            m_setup.get_positive_unstable_dir(),
            Leo::Color(1.0, 0.4, 0.4),
            Leo::Color(1.0, 0.4, 0.4),
            false);

        m_unstable_neg = std::make_unique<UnstableManifoldNeg>(
            std::ref(m_core_ref),
            std::cref(m_setup),
            m_setup.get_negative_unstable_dir(),
            Leo::Color(1.0, 1.0, 0.0),
            Leo::Color(1.0, 1.0, 0.0),
            false);

        m_intersections = std::make_unique<ManifoldIntersections>(
            std::ref(m_core_ref),
            std::cref(*m_unstable_pos),
            std::cref(*m_unstable_neg));
    }

private:
    Lyra::Core3d& m_core_ref;
    const Pcr3bpSetupValues<RMap>& m_setup;

    std::unique_ptr<UnstableManifoldPos> m_unstable_pos;
    std::unique_ptr<UnstableManifoldNeg> m_unstable_neg;

    std::unique_ptr<ManifoldIntersections> m_intersections;
};

class RhezManifolds
{
public:
    RhezManifolds(Lyra::Core3d& core_ref, const Pcr3bpSetupValues<RMap>& setup)
        : m_core_ref(core_ref)
        , m_setup(setup)
        , m_func_ptr_1()
        , m_func_ptr_2()
        , m_func_ptr_3()
        , m_func_ptr_4()
    {}

    void refresh(double h, size_t manifold_select)
    {
        const float thickness = 8e-4f;

        if (manifold_select & 0x1)
        {
            Rhez_U_Alpha_Param_23 param;
            param.setup = m_setup;
            param.h = h;
            param.u_ruler = Leo::Ruler<Real>(-0.999, 0.999, 160, 3);
            param.alpha_ruler = Leo::Ruler<Real>(0, 2*M_PI, 160, 2);
            param.thickness = thickness;
            param.color = Leo::Color(0.4, 0.6, 0.0);

            m_func_ptr_1 = std::make_unique<RenderableUAlpha23>(
                std::ref(m_core_ref), std::cref(param));
        }
        else
        {
            m_func_ptr_1.reset();
        }

        if (manifold_select & 0x2)
        {
            Rhez_U_Alpha_XiEta_23 param;
            param.setup = m_setup;
            param.h = h;
            param.xi_ruler = Leo::Ruler<Real>(-2.5, 2.5, 50, 100);
            param.eta_ruler = Leo::Ruler<Real>(-2.5, 2.5, 50, 100);
            param.thickness = 2*thickness;
            param.color = Leo::Color(0.3, 0.6, 0.9);
            m_func_ptr_2 = std::make_unique<RenderableUXiEta23>(
                std::ref(m_core_ref), std::cref(param));
        }
        else
        {
            m_func_ptr_2.reset();
        }

        if (manifold_select & 0x4)
        {
            Rhez_U_Alpha_Param_23 param;
            param.setup = m_setup;
            param.h = h;
            param.u_ruler = Leo::Ruler<Real>(+1.001, +2.0, 80, 3);
            param.alpha_ruler = Leo::Ruler<Real>(0, 2*M_PI, 160, 2);
            param.thickness = thickness;
            param.color = Leo::Color(0.4, 0.6, 0.0);

            m_func_ptr_3 = std::make_unique<RenderableUAlpha23>(
                std::ref(m_core_ref), std::cref(param));
        }
        else
        {
            m_func_ptr_3.reset();
        }

        if (manifold_select & 0x8)
        {
            Rhez_U_Alpha_Param_23 param;
            param.setup = m_setup;
            param.h = h;
            param.u_ruler = Leo::Ruler<Real>(-1.001, -2.0, 80, 3);
            param.alpha_ruler = Leo::Ruler<Real>(0, 2*M_PI, 160, 2);
            param.thickness = thickness;
            param.color = Leo::Color(0.4, 0.6, 0.0);

            m_func_ptr_4 = std::make_unique<RenderableUAlpha23>(
                std::ref(m_core_ref), std::cref(param));
        }
        else
        {
            m_func_ptr_4.reset();
        }
    }

private:
    Lyra::Core3d& m_core_ref;
    const Pcr3bpSetupValues<RMap>& m_setup;

    std::unique_ptr<RenderableUAlpha23> m_func_ptr_1;
    std::unique_ptr<RenderableUXiEta23> m_func_ptr_2;
    std::unique_ptr<RenderableUAlpha23> m_func_ptr_3;
    std::unique_ptr<RenderableUAlpha23> m_func_ptr_4;
};

class CoreInterior : CoreInteriorBase
{
private:
    Lyra::Core3d& m_core_ref;

    const Pcr3bpSetupValues<RMap> m_setup;
    const Pcr3bpSetupValues<IMap> m_setup_i;

    LineManifolds m_line_manifolds;
    RhezManifolds m_rhez_manifolds;

    using RenderableInterpPositive = RenderableInterp<Pcr3bpRegPoincarePositiveU<RMap>>;
    using RenderableInterpNegative = RenderableInterp<Pcr3bpRegPoincareNegativeU<RMap>>;

    std::unique_ptr<RenderableInterpPositive> m_func_ptr_5;
    std::unique_ptr<RenderableInterpNegative> m_func_ptr_6;
    std::unique_ptr<RenderableInterpPositive> m_func_ptr_7;
    std::unique_ptr<RenderableInterpNegative> m_func_ptr_8;

    RhezUAlpha_Render2<RMap> m_rhez_u_alpha;
    Pcr3bpRegPoincarePositiveU<RMap> m_poincare_positive_u;

    std::unique_ptr<V3> m_v3_1_ptr;
    std::unique_ptr<V3> m_v3_2_ptr;
    std::unique_ptr<V3> m_v3_3_ptr;
    std::unique_ptr<V3> m_v3_4_ptr;
    std::unique_ptr<V3> m_v3_5_ptr;
    std::unique_ptr<V3> m_v3_6_ptr;

public:
    CoreInterior(Lyra::Core3d& core_ref)
        : CoreInteriorBase()
        , m_core_ref(core_ref)
        , m_setup()
        , m_setup_i()

        , m_line_manifolds(m_core_ref, m_setup)
        , m_rhez_manifolds(m_core_ref, m_setup)

        , m_func_ptr_5()
        , m_func_ptr_6()
        , m_func_ptr_7()
        , m_func_ptr_8()

        , m_rhez_u_alpha(m_setup)
        , m_poincare_positive_u(m_setup)
    {
        m_line_manifolds.refresh();
    }

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        CoreInteriorBase::set_param(packet_vector);

        const double h = this->get_param(0);
        const double u_span = this->get_param(1);

        const size_t manifold_select = this->get_param(2);

        const size_t level = this->get_param(3);
        const double points = this->get_param(4);
        const double sub_points = this->get_param(5);
        const double min_multiplier = this->get_param(6);
        const double max_multiplier = this->get_param(7);

        const bool show_interpolation_nodes = this->get_param(8) != 0;

        m_rhez_manifolds.refresh(h, manifold_select);

        {
            InterpParam param;
            param.setup = m_setup;
            param.level = level;
            param.points = points;
            param.sub_points = sub_points;
            param.min_multiplier = exp(min_multiplier)-exp(-30);
            param.max_multiplier = exp(max_multiplier);
            param.line_thickness = 1e-3f;
            param.point_thickness = 1e-2f;
            param.show_interpolation_nodes = show_interpolation_nodes;

            if (manifold_select & 0x10)
            {
                param.direction = m_setup.get_positive_unstable_dir();
                param.color = Leo::Color(1.0, 0.0, 0.0);

                m_func_ptr_5 = std::make_unique<RenderableInterpPositive>(std::ref(m_core_ref), std::cref(param));
            }
            else
            {
                m_func_ptr_5.reset();
            }

            if (manifold_select & 0x20)
            {
                param.direction = m_setup.get_negative_unstable_dir();
                param.color = Leo::Color(1.0, 1.0, 0.0);

                m_func_ptr_6 = std::make_unique<RenderableInterpNegative>(std::ref(m_core_ref), std::cref(param));
            }
            else
            {
                m_func_ptr_6.reset();
            }

            if (manifold_select & 0x40)
            {
                param.direction = -m_setup.get_positive_unstable_dir();
                param.color = Leo::Color(1.0, 0.0, 0.0);

                m_func_ptr_7 = std::make_unique<RenderableInterpPositive>(std::ref(m_core_ref), std::cref(param));
            }
            else
            {
                m_func_ptr_7.reset();
            }

            if (manifold_select & 0x80)
            {
                param.direction = -m_setup.get_negative_unstable_dir();
                param.color = Leo::Color(1.0, 1.0, 0.0);

                m_func_ptr_8 = std::make_unique<RenderableInterpNegative>(std::ref(m_core_ref), std::cref(param));
            }
            else
            {
                m_func_ptr_8.reset();
            }
        }

        {
            const double x1 = this->get_param(9);
            const double y1 = this->get_param(10);
            const double z1 = this->get_param(11);

            const RVector dir = m_setup.get_positive_unstable_dir();

            RMatrix dummy1;
            const RVector e = m_rhez_u_alpha(dir * x1, dummy1);

            RMatrix dummy2;
            const RVector g = m_poincare_positive_u(e);
            std::cout << g << '\n';

            const RVector h = m_poincare_positive_u(g);
            const RVector i = m_poincare_positive_u(h);
            const RVector j = m_poincare_positive_u(i);
            const RVector k = m_poincare_positive_u(j);

            m_v3_1_ptr = std::make_unique<V3>(std::ref(m_core_ref), RVector{ e[2], e[3], e[0] });
            m_v3_2_ptr = std::make_unique<V3>(std::ref(m_core_ref), RVector{ g[2], g[3], g[0] });
            m_v3_3_ptr = std::make_unique<V3>(std::ref(m_core_ref), RVector{ h[2], h[3], h[0] });
            m_v3_4_ptr = std::make_unique<V3>(std::ref(m_core_ref), RVector{ i[2], i[3], i[0] });
            m_v3_5_ptr = std::make_unique<V3>(std::ref(m_core_ref), RVector{ j[2], j[3], j[0] });
            m_v3_6_ptr = std::make_unique<V3>(std::ref(m_core_ref), RVector{ k[2], k[3], k[0] });
        }
    }
};

}

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore3d<Ursa::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Ursa::create_window_properties("RHEZ U"), argc, argv);

    window.show();
    window.run();
    return 0;
}
