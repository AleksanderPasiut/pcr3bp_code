///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "tools/test_tools.hpp"
#include <carina/enp_map.hpp>
#include "tools/interpolation_nodes_renderable.hpp"
#include "capd_renderable.hpp"

#include <leo/interpolation/interpolation.hpp>

#include "ppv_manifolds.hpp"

namespace Ursa
{

struct InterpParam
{
    Pcr3bpSetupValues<RMap> setup;

    RVector direction;

    size_t level;
    size_t points;
    size_t sub_points;
    double min_multiplier;
    double max_multiplier;
    Leo::Color color;

    float line_thickness;
    float point_thickness;

    bool show_interpolation_nodes;
};

template<typename PoincareT>
class RenderableInterp
{
public:
    RenderableInterp(Lyra::Core3d& core_ref, const InterpParam& param)
        : m_core_ref(core_ref)
        , m_factory(param.setup)
        , m_ppv()
        , m_ppv_proj()
        , m_renderables_line()
    {
        PpvManifoldInitParam fparam;
        fparam.origin = { 0.0, 0.0 };
        fparam.direction = param.direction;
        fparam.points = param.points;
        fparam.min_multiplier = param.min_multiplier;
        fparam.max_multiplier = param.max_multiplier;
        m_ppv.emplace_back( m_factory.init(fparam) );
        m_ppv_proj.emplace_back( std::ref(*m_ppv.rbegin()) );

        for (size_t i = 0; i < param.level; ++i)
        {
            PpvManifold<RMap>& ref = *m_ppv.rbegin();
            m_ppv.emplace_back( m_factory.evolve(ref) );
            m_ppv_proj.emplace_back( std::ref(*m_ppv.rbegin()) );
        }

        init_renderables(param);
    }

    virtual ~RenderableInterp() noexcept
    {
        for (auto & renderable : m_renderables_line)
        {
            m_core_ref.unregister_manifold(&renderable);
        }

        for (auto & renderable : m_renderables_points)
        {
            m_core_ref.unregister_manifold(renderable.get_renderable_ptr());
        }
    }

private:
    void init_renderables(const InterpParam& param)
    {
        const Leo::RulerSet<1> ruler_set = { Leo::Ruler<>(
            m_ppv.rbegin()->get_min_argument(),
            m_ppv.rbegin()->get_max_argument(),
            m_ppv.rbegin()->get_nodes().size(),
            param.sub_points ) };

        for (auto & ppv : m_ppv_proj)
        {
            m_renderables_line.emplace_back( std::ref(ppv), ruler_set, param.color);
        }

        for (auto & renderable : m_renderables_line)
        {
            renderable.fill(param.line_thickness);
            m_core_ref.register_manifold(&renderable);
        }

        if (param.show_interpolation_nodes)
        {
            for (auto & ppv : m_ppv)
            {
                m_renderables_points.emplace_back(
                    std::cref(ppv.get_nodes()),
                    5,
                    Carina::IdxList<size_t>{ 2, 3, 0 },
                    param.color,
                    param.point_thickness);
            }

            for (auto & renderable : m_renderables_points)
            {
                m_core_ref.register_manifold(renderable.get_renderable_ptr());
            }
        }
    }

    Lyra::Core3d& m_core_ref;

    PpvManifoldFactory<PoincareT> m_factory;

    class PpvManifoldWithProjection
    {
    public:
        using VectorType = RVector;

        PpvManifoldWithProjection(PpvManifold<RMap>& ppv)
            : m_ppv(ppv)
            , m_projection(Carina::ProjectionMap<RMap>::create(5, {2, 3, 0}))
        {}

        RVector operator() (const RVector& arg)
        {
            return m_projection(m_ppv(arg));
        }

    private:
        PpvManifold<RMap>& m_ppv;
        RMap m_projection;
    };

public:
    std::list<PpvManifold<RMap>> m_ppv;
    std::list<PpvManifoldWithProjection> m_ppv_proj;
    std::list<InterpolationNodesRenderable<RMap, 3>> m_renderables_points;
    std::list<CapdMapRenderable<PpvManifoldWithProjection, 1, 3>> m_renderables_line;
};

}
