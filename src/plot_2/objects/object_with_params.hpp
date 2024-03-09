///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <memory>

namespace Pcr3bpProof
{

template<typename ObjectT>
class Renderable4d_WithParams
{
public:
    using ParamsT = typename ObjectT::Params;

    Renderable4d_WithParams(
        Lyra::Core3dObjects& core_objects_ref,
        const Manifold4_Transformation& transformation_ref)
            : m_core_objects_ref(core_objects_ref)
            , m_transformation_ref(transformation_ref)
    {}

    void show(const ParamsT& params)
    {
        if (m_object_ptr)
        {
            if (m_object_ptr->get_params() != params)
            {
                rebuild_object(params);
            }
        }
        else
        {
            rebuild_object(params);
        }
    }

    void hide()
    {
        m_object_ptr.reset();
    }

    void refresh()
    {
        if (m_object_ptr)
        {
            m_object_ptr->refresh();
        }
    }

private:
    void rebuild_object(const ParamsT& params)
    {
        m_object_ptr = std::make_unique<ObjectT>(
            std::ref(m_core_objects_ref),
            std::cref(m_transformation_ref),
            std::cref(params));
    }

    Lyra::Core3dObjects& m_core_objects_ref;
    const Manifold4_Transformation& m_transformation_ref;

    std::unique_ptr<ObjectT> m_object_ptr {};
};

}
