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
private:
    using Vector4 = CapdVectorRenderable4<RVector>;

    Lyra::Core3d& m_core_ref;

    Vector4 m_vector;

public:
    HL_Map(Lyra::Core3d& core_ref, RVector U, const Leo::Matrix4f & matrix)
        : m_core_ref(core_ref)
        , m_vector(core_ref.get_objects(), RVector{ U[2], U[3], U[0], U[1] }, std::cref(matrix))
    {
        refresh();

        // m_core_ref.register_manifold(&m_vector);
    }

    virtual ~HL_Map() noexcept
    {
        // m_core_ref.unregister_manifold(&m_vector);
    }

    void refresh()
    {
        m_vector.fill(0.05f);
    }
};

}
