///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_2d.hpp>
#include <lyra/core2d.hpp>

#include "plot_common/window_properties.hpp"
#include "plot_common/core_interior_base.hpp"

#include "pcr3bp_obsolete/pcr3bp_reg_lyapunov_orbit_implicit.hpp"
#include "pcr3bp_obsolete/pcr3bp_reg_poincare.hpp"

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

        const Real u0 = this->get_param(0);
        const size_t point_count = static_cast<size_t>(this->get_param(1));
        const size_t steps = static_cast<size_t>(this->get_param(2));

        const RVector PV = LyapunovOrbitRegParam::calculate(u0, steps);

        RegEvolutionParam param;
        param.setup = m_masses.get_setup();
        param.u0 = PV[0];
        param.v0 = PV[1];
        param.pu0 = PV[2];
        param.pv0 = PV[3];
        param.h = PV[4];

        Pcr3bpRegPoincarePositiveU<RMap> poincare(m_masses.get_setup());
        Real t = poincare.get_return_time(RVector{ param.u0, param.v0, param.pu0, param.pv0, param.h });
        param.t = 2 * t;

        param.point_count = point_count;

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
        Ursa::create_window_properties("reg lyapunov orbit"), argc, argv);

    window.show();
    window.run();
    return 0;
}
