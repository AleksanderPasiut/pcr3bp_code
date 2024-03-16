///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "capd_renderable.hpp"
#include "tools/test_tools.hpp"

namespace Pcr3bpProof
{

class HL_Map
{
public:
    struct Params
    {
        RVector U;
        float thickness;

        bool operator!= (const Params& arg)
        {
            return 
                arg.thickness != this->thickness ||
                arg.U != this->U;
        }
    };

    HL_Map(
        Lyra::Core3dObjects& core_objects_ref,
        const Manifold4_Transformation & transformation_ref,
        const Params& params)
            : m_params(params)
            , m_vector(core_objects_ref, m_params.U, std::cref(transformation_ref))
    {
        refresh();
    }

    void refresh()
    {
        m_vector.fill(m_params.thickness);
    }

private:
    using Vector4 = CapdVectorRenderable4<RVector>;

    Params m_params;

    Vector4 m_vector;
};

}
