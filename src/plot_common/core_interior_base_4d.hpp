///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <leo/matrix4f.hpp>
#include <leo/array/array_access_std.hpp>
#include <lyra/manifold/manifold4.hpp>
#include "plot_common/core_interior_base.hpp"

#include <iostream>

namespace Pcr3bpProof
{

class CoreInteriorBase4d : public CoreInteriorBase
{
public:
    CoreInteriorBase4d() : CoreInteriorBase()
    {
        m_rotation.set_identity();

        m_transformation = [this](Lyra::Point4d arg) -> Lyra::Point4d
        {
            Leo::ArrayAccessStd<4, float> in(arg);

            Leo::ArrayAccessStdRead offset_access(m_offset);
            Leo::ArrayElemOp::sub(in, offset_access);
            Leo::ArrayScalarOp::mul(in, m_scale);

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

    void set_offset(Lyra::Point4d offset)
    {
        m_offset = { offset[2], offset[3], offset[0], offset[1] };
    }

    Lyra::Manifold4_Transformation const & get_transformation() const noexcept
    {
        return m_transformation;
    }

private:
    Leo::Matrix4f m_rotation {};
    float m_scale { 1.0f };
    Lyra::Point4d m_offset {};

    Lyra::Manifold4_Transformation m_transformation {};
};

}
