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
        Leo::ArrayFill::mov(m_reorder, 
        {
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0
        });

        Leo::ArrayFill::mov(m_rotation, m_reorder);
    }

    void set_rotation(Leo::Matrix4f const & matrix) noexcept
    {
        Leo::ArrayTensorOp::mul<1, 0>(m_rotation, matrix, m_reorder);
        m_transformation.inc_params_counter();
    }

    void set_scale(double scale) noexcept
    {
        m_scale = scale;
        m_transformation.inc_params_counter();
    }

    void set_offset(std::array<double, 4> offset)
    {
        m_offset = offset;
        m_transformation.inc_params_counter();
    }

    Manifold4_Transformation const & get_transformation() const noexcept
    {
        return m_transformation;
    }

private:
    Leo::Matrix4f m_reorder {};

    Leo::Matrix4f m_rotation {};
    double m_scale { 1.0f };
    std::array<double, 4> m_offset {};

    Manifold4_Transformation m_transformation
    {
        [this](std::array<double, 4> arg) -> Lyra::Point4d
        {
            Leo::ArrayAccessStd<4, double> in(arg);

            Leo::ArrayAccessStdRead offset_access(m_offset);
            Leo::ArrayElemOp::sub(in, offset_access);
            Leo::ArrayScalarOp::mul(in, m_scale);

            Lyra::Point4d ret {};
            Leo::ArrayAccessStd out(ret);
            Leo::ArrayTensorOp::mul<1, 0>(out, m_rotation, in);
            return ret;
        }
    };
};

}
