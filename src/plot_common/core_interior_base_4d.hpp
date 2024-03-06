///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <leo/matrix4f.hpp>
#include <leo/array/array_access_std.hpp>
#include <lyra/manifold/manifold4.hpp>
#include "plot_common/core_interior_base.hpp"


namespace Pcr3bpProof
{

class CoreInteriorBase4d : public CoreInteriorBase
{
public:
    CoreInteriorBase4d() : CoreInteriorBase(), m_rotation()
    {
        m_rotation.set_identity();

        m_transformation = [this](Lyra::Point4d const & arg) -> Lyra::Point4d
        {
            Leo::ArrayAccessStdRead<4, float> in(arg);

            Lyra::Point4d ret {};
            Leo::ArrayAccessStd out(ret);

            Leo::ArrayTensorOp::mul<1, 0>(out, m_rotation, in);
            return ret;
        };
    }

    void set_rotation(Leo::Matrix4f const & matrix) noexcept
    {
        m_rotation = matrix;
    }

    void set_scale(float scale) noexcept
    {
        m_scale = scale;
    }

    Lyra::Manifold4_Transformation const & get_transformation() const noexcept
    {
        return m_transformation;
    }

private:
    Leo::Matrix4f m_rotation;
    float m_scale { 1.0f };

    Lyra::Manifold4_Transformation m_transformation {};
};

}
