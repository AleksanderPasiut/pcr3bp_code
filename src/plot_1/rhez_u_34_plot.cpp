///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_4d.hpp>
#include <lyra/core3d.hpp>

#include "plot_common/window_properties.hpp"
#include "plot_common/core_interior_base_4d.hpp"

#include "capd_renderable.hpp"

#include "plot_1/objects/rhez_u_alpha_object_34.hpp"

namespace Ursa
{

class CoreInterior : CoreInteriorBase4d
{
private:
    Lyra::Core3d& m_core_ref;

    RenderableUAlpha34 m_renderable;

public:
    CoreInterior(Lyra::Core3d& core_ref)
        : CoreInteriorBase4d()
        , m_core_ref(core_ref)
        , m_renderable(m_core_ref, get_param(), this->get_rotation())
    {}

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        CoreInteriorBase::set_param(packet_vector);
    }

    void set_rotation_4d(Leo::Matrix4f const & matrix)
    {
        CoreInteriorBase4d::set_rotation(matrix);
        m_renderable.refresh();
    }

    static Rhez_U_Alpha_Param_34 get_param()
    {
        Rhez_U_Alpha_Param_34 param;
        param.setup = Pcr3bpSetupValues<RMap>();
        param.h_ruler = Leo::Ruler<double>(-1.29, +1.29, 20, 10);
        param.u_ruler = Leo::Ruler<double>(-0.997, +0.997, 10, 10);
        param.alpha_ruler = Leo::Ruler<double>(0.0, 2*M_PI, 20, 10);
        param.thickness = 3e-3;
        return param;
    }
};

}

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore4d<Ursa::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Ursa::create_window_properties("RHEZ U 34"), argc, argv);

    window.show();
    window.run();
    return 0;
}
