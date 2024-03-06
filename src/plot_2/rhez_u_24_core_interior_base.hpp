///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "plot_common/core_interior_base_4d.hpp"

#include "capd_renderable.hpp"

namespace Pcr3bpProof
{

class CoreInteriorBaseRhez_u_24 : protected CoreInteriorBase4d
{
private:
    Lyra::Core3d& m_core_ref;

protected:
    CoreInteriorBaseRhez_u_24(Lyra::Core3d& core_ref)
        : CoreInteriorBase4d()
        , m_core_ref(core_ref)
    {}

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        CoreInteriorBase4d::set_param(packet_vector);
    }

    void set_rotation_4d(Leo::Matrix4f const & matrix)
    {
        CoreInteriorBase4d::set_rotation(matrix);
    }

    Lyra::Core3d& get_core_ref()
    {
        return m_core_ref;
    }
};

}
