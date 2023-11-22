///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <leo/matrix4f.hpp>
#include "plot_common/core_interior_base.hpp"

namespace Ursa
{

class CoreInteriorBase4d : public CoreInteriorBase
{
public:
    CoreInteriorBase4d() : CoreInteriorBase(), m_rotation()
    {
        m_rotation.set_identity();
    }

    void set_rotation(Leo::Matrix4f const & matrix) noexcept
    {
        m_rotation = matrix;
    }

    const Leo::Matrix4f& get_rotation() const noexcept
    {
        return m_rotation;
    }

private:
    Leo::Matrix4f m_rotation;
};

}
