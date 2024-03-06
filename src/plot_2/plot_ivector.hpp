///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>
#include <lyra/manifold/manifold_interface.hpp>
#include "tools/test_tools.hpp"

namespace Pcr3bpProof
{

class PlotIVector
{
public:
    struct Params
    {
        float thickness {};
        Leo::Color color {};

        double x1_scale {};
        double x2_scale {};
    };

    PlotIVector(Lyra::Core2d& core_ref, const IVector& iv, const Params& params)
        : m_core_ref(core_ref)
        , m_renderable({ Leo::RulerSet<2>({
            Leo::Ruler<double>(iv[0].leftBound(), iv[0].rightBound(), 2, 1),
            Leo::Ruler<double>(iv[1].leftBound(), iv[1].rightBound(), 2, 1),
        }), params.color })
    {
        auto func = [this, params](const std::array<double, 2>& in) -> std::array<float, 2>
        {
            return {
                float( params.x1_scale * in[0] ),
                float( params.x2_scale * in[1] )
            };
        };

        m_renderable.fill( func, params.thickness );
        m_core_ref.register_manifold(&m_renderable);
    }

    virtual ~PlotIVector() noexcept
    {
        m_core_ref.unregister_manifold(&m_renderable);
    }

private:
    Lyra::Core2d& m_core_ref;
    Lyra::ManifoldInterface<2,2> m_renderable;
};

}
