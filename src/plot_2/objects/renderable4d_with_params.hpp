///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "manifold4_transformation.hpp"
#include <memory>

#include "plot_common/timelevel_divisor.hpp"
#include "plot_common/edge_detector.hpp"

namespace Pcr3bpProof
{

template<typename ObjectT>
class Renderable4d_WithParams
{
public:
    using ParamsT = typename ObjectT::Params;

    Renderable4d_WithParams(
        Lyra::Core3d& core_ref,
        const Manifold4_Transformation& transformation_ref)
            : m_core_objects_ref(core_ref.get_objects())
            , m_transformation_ref(transformation_ref)
    {}

    void rebuild(const ParamsT& params)
    {
        if (m_object_ptr)
        {
            if (m_object_ptr->get_params() != params)
            {
                rebuild_object(params);
                m_divisor.reset();
                m_visible = true;
                m_edge_detector.reset(true);
            }
            else if (m_transformation_ref.update())
            {
                m_object_ptr->refresh();
            }
        }
        else
        {
            rebuild_object(params);
            m_divisor.reset();
            m_visible = true;
            m_edge_detector.reset(true);
        }
    }

    void hide()
    {
        m_object_ptr.reset();
        m_divisor.reset();
        m_visible = false;
        m_edge_detector.reset(false);
    }

    void highlight(bool arg)
    {
        m_highlighted = arg;
    }

    void refresh()
    {
        if (m_object_ptr)
        {
            if (m_transformation_ref.update())
            {
                m_object_ptr->refresh();
            }
        }
    }

    void heartbeat()
    {
        if (m_highlighted)
        {
            bool const div_output = m_divisor.update();

            if (div_output)
            {
                m_visible = !m_visible;
            }
        }
        else
        {
            m_visible = true;
        }

        if (m_edge_detector.update(m_visible))
        {
            if (m_object_ptr)
            {
                m_object_ptr->show(m_visible);
            }
        }
    }

private:
    void rebuild_object(const ParamsT& params)
    {
        m_object_ptr = std::make_unique<ObjectT>(
            std::ref(m_core_objects_ref),
            std::cref(m_transformation_ref.ref()),
            std::cref(params));
    }

    Lyra::Core3dObjects& m_core_objects_ref;
    Manifold4_Transformation_Ref m_transformation_ref;

    std::unique_ptr<ObjectT> m_object_ptr {};

    TimelevelDivisor m_divisor {};

    bool m_visible {};
    bool m_highlighted {};

    EdgeDetector m_edge_detector {};
};

}
