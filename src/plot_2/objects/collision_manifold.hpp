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
    struct Param
    {
        float thickness;
    };

    CollisionManifold(
        Lyra::Core3dObjects& core_objects_ref,
        const Manifold4_Transformation & transformation_ref,
        const Param& param)
            : m_renderable(
                core_objects_ref,
                m_map,
                Leo::RulerSet<1>{ Leo::Ruler( 0.0, 2 * M_PI, 10000, 1 ) },
                std::cref(transformation_ref))
            , m_param(param)
    {
        refresh();
    }

    void refresh()
    {
        m_renderable.fill(m_param.thickness);
    }

private:
    static RMap create_collision_manifold_map()
    {
        using CapdUtils::Node;
        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& phi = in[0];
            Node& r = param[0];

            out[2] = Node(0.0);
            out[3] = Node(0.0);
            out[0] = r * sin(phi);
            out[1] = r * cos(phi);
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

    Param m_param;
};

}
