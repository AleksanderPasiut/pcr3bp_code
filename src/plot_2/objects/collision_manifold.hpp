///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "capd_renderable.hpp"
#include "tools/test_tools.hpp"

namespace Pcr3bpProof
{

class CollisionManifold
{
public:
    struct Params
    {
        float thickness;

        bool operator!= (const Params& params) const noexcept
        {
            return params.thickness != thickness;
        }
    };

    CollisionManifold(
        Lyra::Core3dObjects& core_objects_ref,
        const Manifold4_Transformation & transformation_ref,
        const Params& params)
            : m_renderable(
                core_objects_ref,
                m_map,
                Leo::RulerSet<1>{ Leo::Ruler( 0.0, 2 * M_PI, 10000, 1 ) },
                std::cref(transformation_ref))
            , m_params(params)
    {
        refresh();
    }

    void refresh()
    {
        m_renderable.fill(m_params.thickness);
    }

    const Params& get_params()
    {
        return m_params;
    }

    void show(bool arg)
    {
        m_renderable.show(arg);
    }

private:
    static RMap create_collision_manifold_map()
    {
        using CapdUtils::Node;
        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& phi = in[0];
            Node& r = param[0];

            out[0] = Node(0.0);
            out[1] = Node(0.0);
            out[2] = r * sin(phi);
            out[3] = r * cos(phi);
        };

        RMap map(func, 1, 4, 1);
        map.setParameter(0, 2.81112771399491);
        return map;
    }

    RMap m_map
    {
        create_collision_manifold_map()
    };

    CapdMapRenderable<decltype(m_map), 1, 4> m_renderable;

    Params m_params;
};

}
