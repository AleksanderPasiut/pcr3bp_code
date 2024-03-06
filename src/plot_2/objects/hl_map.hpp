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
    HL_Map(Lyra::Core3d& core_ref, RVector U, const Leo::Matrix4f & matrix)
        : m_core_ref(core_ref)
        , m_vector(core_ref.get_objects(), RVector{ U[2], U[3], U[0], U[1] }, std::cref(matrix))
    {
        refresh();
    }

    void refresh()
    {
        m_vector.fill(0.05f);
    }

private:
    using Vector4 = CapdVectorRenderable4<RVector>;

    Lyra::Core3d& m_core_ref;

    Vector4 m_vector;
};

}
