///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>

#include "capd_renderable.hpp"
#include "tools/test_tools.hpp"

namespace Pcr3bpProof
{
    template<typename MapT>
    class Plot22 final
    {
    private:
        Lyra::Core2d& m_core_ref;
        MapT & m_map;
        Lyra::ManifoldInterface<2,2> m_renderable;

        using Map22v = CapdVectorRenderable<CapdUtils::RVector, 2>;

        std::list<Map22v> m_map22v;

        Lyra::ManifoldInterface<1,2> m_renderable_1;
        Lyra::ManifoldInterface<1,2> m_renderable_2;

        const double m_x1_scale;
        const double m_x2_scale;

    public:
        struct Params
        {
            Leo::Ruler<Real> x1_ruler;
            Leo::Ruler<Real> x2_ruler;
            float thickness;

            double x1_scale;
            double x2_scale;
        };

        Plot22(Lyra::Core2d& core_ref, MapT& map, const Params& params)
            : m_core_ref(core_ref)
            , m_map(map)
            , m_renderable({ Leo::RulerSet<2>({ params.x1_ruler, params.x2_ruler }), Leo::Color(0.7, 0.8, 0.3) })
            , m_map22v()
            , m_renderable_1({ Leo::RulerSet<1>({ params.x1_ruler }), Leo::Color(1.0, 0.0, 0.3) })
            , m_renderable_2({ Leo::RulerSet<1>({ params.x2_ruler }), Leo::Color(0.0, 0.3, 1.0) })
            , m_x1_scale(params.x1_scale)
            , m_x2_scale(params.x2_scale)
        {
            fill_renderable(params.thickness);

            const RVector x0 =
            {
                0.5 * (params.x1_ruler.maximum() + params.x1_ruler.minimum()),
                0.5 * (params.x2_ruler.maximum() + params.x2_ruler.minimum())
            };

            fill_renderable_1(x0, params.thickness);
            fill_renderable_2(x0, params.thickness);
            fill_map22v(x0, params);
        }

        virtual ~Plot22() noexcept
        {
            for (Map22v & m : m_map22v)
            {
                m_core_ref.unregister_manifold(&m);
            }

            m_map22v.clear();

            m_core_ref.unregister_manifold(&m_renderable);

            m_core_ref.unregister_manifold(&m_renderable_1);
            m_core_ref.unregister_manifold(&m_renderable_2);
        }

    private:
        RVector evaluate_internal(RVector x)
        {
            RVector ret = m_map(x);
            return { m_x1_scale * ret[0], m_x2_scale * ret[1] };
        }

        std::array<float, 2> evaluate(double x1, double x2)
        {
            const RVector X = { x1, x2 };
            const RVector Y = evaluate_internal(X);
            return { float(Y[0]), float(Y[1]) };
        }

        void fill_renderable(float thickness)
        {
            auto func = [this](const std::array<double, 2>& in) -> std::array<float, 2>
            {
                return evaluate( in[0], in[1] );
            };

            m_renderable.fill( func, thickness );
            m_core_ref.register_manifold(&m_renderable);
        }

        void fill_renderable_1(const RVector x0, float thickness)
        {
            auto func = [this, x0](const std::array<double, 1>& in) -> std::array<float, 2>
            {
                return evaluate( in[0], x0[1] );
            };

            m_renderable_1.fill( func, thickness );
            m_core_ref.register_manifold(&m_renderable_1);
        }

        void fill_renderable_2(const RVector x0, float thickness)
        {
            auto func = [this, x0](const std::array<double, 1>& in) -> std::array<float, 2>
            {
                return evaluate( x0[0], in[0] );
            };

            m_renderable_2.fill( func, thickness );
            m_core_ref.register_manifold(&m_renderable_2);
        }

        void fill_map22v(const RVector x0, const Params& params)
        {
            const RVector x1 = { params.x1_ruler.maximum(), x0[1] };
            const RVector x2 = { x0[0], params.x2_ruler.maximum() };

            m_map22v.emplace_back( evaluate_internal(x0), Leo::Color(0.3, 1.0, 0.0) );
            m_map22v.emplace_back( evaluate_internal(x1), Leo::Color(1.0, 0.0, 0.3) );
            m_map22v.emplace_back( evaluate_internal(x2), Leo::Color(0.0, 0.3, 1.0) );

            for (Map22v & m : m_map22v)
            {
                m.fill(params.thickness * 5);
                m_core_ref.register_manifold(&m);
            }
        }
    };
}
