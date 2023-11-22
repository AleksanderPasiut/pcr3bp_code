///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"
#include "capd_renderable.hpp"

namespace Ursa
{

class V3
{
public:
    V3(Lyra::Core3d& core_ref, const RVector& v)
        : m_core_ref(core_ref)
        , m_vector(v, Leo::Color(0.0, 1.0, 0.0))
    {
        m_vector.fill(1e-2f);
        m_core_ref.register_manifold(&m_vector);
    }

    ~V3() noexcept
    {
        m_core_ref.unregister_manifold(&m_vector);
    }

private:
    Lyra::Core3d& m_core_ref;

    CapdVectorRenderable<RVector, 3> m_vector;
};

}
