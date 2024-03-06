///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "plot_common/core_interior_base_4d.hpp"

#include "capd_renderable.hpp"

#if 0

#include "plot_1/objects/rhez_u_alpha_object_24.hpp"

#endif

namespace Pcr3bpProof
{

class CoreInteriorBaseRhez_u_24 : protected CoreInteriorBase4d
{
private:
    Lyra::Core3d& m_core_ref;

#if 0
    std::unique_ptr<RenderableUAlpha24> m_manifold_mid;
    std::unique_ptr<RenderableUAlpha24> m_manifold_pos;
    std::unique_ptr<RenderableUAlpha24> m_manifold_neg;
#endif

protected:
    CoreInteriorBaseRhez_u_24(Lyra::Core3d& core_ref)
        : CoreInteriorBase4d()
        , m_core_ref(core_ref)
    {}

#if 0

    void reload_pos_manifold(double h, bool show, float thickness)
    {
        if (show)
        {
            Rhez_U_Alpha_Param_24 param;
            param.setup = Pcr3bpSetupValues<RMap>();
            param.h = h;
            param.u_ruler = Leo::Ruler<Real>(1.001, 2.0, 20, 2);
            param.alpha_ruler = Leo::Ruler<Real>(0, 2*M_PI, 60, 2);
            param.thickness = thickness;

            m_manifold_pos = std::make_unique<RenderableUAlpha24>(
                std::ref(m_core_ref), std::cref(param), std::cref(this->get_rotation()));
        }
        else
        {
            m_manifold_pos.reset();
        }
    }

    void reload_mid_manifold(double h, bool show, float thickness)
    {
        if (show)
        {
            Rhez_U_Alpha_Param_24 param;
            param.setup = Pcr3bpSetupValues<RMap>();
            param.h = h;
            param.u_ruler = Leo::Ruler<Real>(-0.999, 0.999, 40, 2);
            param.alpha_ruler = Leo::Ruler<Real>(0, 2*M_PI, 60, 2);
            param.thickness = thickness;

            m_manifold_mid = std::make_unique<RenderableUAlpha24>(
                std::ref(m_core_ref), std::cref(param), std::cref(this->get_rotation()));
        }
        else
        {
            m_manifold_mid.reset();
        }
    }

    void reload_neg_manifold(double h, bool show, float thickness)
    {
        if (show)
        {
            Rhez_U_Alpha_Param_24 param;
            param.setup = Pcr3bpSetupValues<RMap>();
            param.h = h;
            param.u_ruler = Leo::Ruler<Real>(-1.001, -2.0, 20, 2);
            param.alpha_ruler = Leo::Ruler<Real>(0, 2*M_PI, 60, 2);
            param.thickness = thickness;

            m_manifold_neg = std::make_unique<RenderableUAlpha24>(
                std::ref(m_core_ref), std::cref(param), std::cref(this->get_rotation()));
        }
        else
        {
            m_manifold_neg.reset();
        }
    }

#endif

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        CoreInteriorBase4d::set_param(packet_vector);
    }

    void set_rotation_4d(Leo::Matrix4f const & matrix)
    {
        CoreInteriorBase4d::set_rotation(matrix);
    
#if 0

        if (m_manifold_mid)
        {
            m_manifold_mid->refresh();
        }

        if (m_manifold_pos)
        {
            m_manifold_pos->refresh();
        }

        if (m_manifold_neg)
        {
            m_manifold_neg->refresh();
        }

#endif
    }

    Lyra::Core3d& get_core_ref()
    {
        return m_core_ref;
    }
};

}
