///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_2d.hpp>
#include <lyra/core2d.hpp>

#include "plot_common/window_properties.hpp"
#include "plot_common/core_interior_base.hpp"

#include "objects/reg_masses.hpp"
#include "objects/reg_evolution.hpp"

namespace Ursa
{

class CoreInterior : CoreInteriorBase
{
private:
    Lyra::Core2d& m_core_ref;

    RegMasses m_masses;

    std::unique_ptr<RegEvolution> m_func_ptr;

public:
    CoreInterior(Lyra::Core2d& core_ref)
        : CoreInteriorBase()
        , m_core_ref(core_ref)
        , m_masses(core_ref)
        , m_func_ptr()
    {}

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        CoreInteriorBase::set_param(packet_vector);

        RegEvolutionParam param;
        param.setup = m_masses.get_setup();
        param.u0 = this->get_param(0);
        param.v0 = this->get_param(1);
        param.pu0 = this->get_param(2);
        param.pv0 = this->get_param(3);
        param.h = this->get_param(4);
        param.t = this->get_param(5);
        param.point_count = this->get_param(6);

        param.point_thickness = 5e-3f;
        param.line_thickness = 1e-3f;
        param.point_subcount = 10;
        param.color = Leo::Color(0.3, 0.1, 0.8);

        m_func_ptr = std::make_unique<RegEvolution>(std::ref(m_core_ref), std::cref(param));
    }
};

}

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore2d<Ursa::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Ursa::create_window_properties("reg evolution raw"), argc, argv);

    window.show();
    window.run();
    return 0;
}
