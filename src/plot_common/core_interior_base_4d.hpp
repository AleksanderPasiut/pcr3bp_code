///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <leo/matrix4f.hpp>
#include <leo/array/array_access_std.hpp>
#include <lyra/manifold/manifold4.hpp>
#include "plot_common/core_interior_base.hpp"

#include "capd_renderable.hpp"

#include <iostream>

namespace Pcr3bpProof
{

class CoreInteriorBase4d : public CoreInteriorBase
{
public:
    CoreInteriorBase4d() : CoreInteriorBase()
    {
        m_rotation.set_identity();

        m_transformation = [this](std::array<double, 4> arg) -> Lyra::Point4d
        {
            Leo::ArrayAccessStd<4, double> in(arg);

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

    void set_scale(double scale) noexcept
    {
        m_scale = scale;
    }

    void set_offset(std::array<double, 4> offset)
    {
        m_offset = { offset[2], offset[3], offset[0], offset[1] };
    }

    Manifold4_Transformation const & get_transformation() const noexcept
    {
        return m_transformation;
    }

private:
    Leo::Matrix4f m_rotation {};
    double m_scale { 1.0f };
    std::array<double, 4> m_offset {};

    Manifold4_Transformation m_transformation {};
};

}
