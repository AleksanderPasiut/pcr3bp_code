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
    struct Param
    {
        RVector U;
        float thickness;

        bool operator!= (const Param& arg)
        {
            return 
                arg.thickness != this->thickness ||
                arg.U != this->U;
        }
    };

    HL_Map(
        Lyra::Core3dObjects& core_objects_ref,
        const Manifold4_Transformation & transformation_ref,
        const Param& param)
            : m_param(param)
            , m_vector(
                core_objects_ref,
                RVector{ m_param.U[2], m_param.U[3], m_param.U[0], m_param.U[1] },
                std::cref(transformation_ref))
    {
        refresh();
    }

    void refresh()
    {
        m_vector.fill(m_param.thickness);
    }

private:
    using Vector4 = CapdVectorRenderable4<RVector>;

    Param m_param;

    Vector4 m_vector;
};

}
