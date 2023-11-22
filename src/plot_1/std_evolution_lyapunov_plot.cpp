///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_2d.hpp>
#include <lyra/core2d.hpp>

#include "plot_common/window_properties.hpp"
#include "plot_common/core_interior_base.hpp"

#include "pcr3bp_obsolete/pcr3bp_std_lyapunov_orbit_implicit.hpp"

#include "objects/std_masses.hpp"
#include "objects/std_evolution.hpp"

namespace Ursa
{

class CoreInterior : CoreInteriorBase
{
private:
    Lyra::Core2d& m_core_ref;

    StdMasses m_masses;

    std::unique_ptr<StdEvolution> m_func_ptr;

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

        const Real x0 = this->get_param(0);
        const size_t point_count = static_cast<size_t>(this->get_param(1));
        const size_t steps = static_cast<size_t>(this->get_param(2));

        Real t {};
        const IVector PY = Pcr3bpStdLyapunovOrbitImplicit::calculate_py(x0, steps, t);
        t *= 2;

        StdEvolutionParam param;
        param.setup = m_masses.get_setup();
        param.x0 = x0;
        param.y0 = 0.0;
        param.px0 = 0.0;
        param.py0 = PY[0].mid().leftBound();
        param.t = t;
        param.point_count = point_count;

        param.point_thickness = 5e-3f;
        param.line_thickness = 1e-3f;
        param.point_subcount = 10;
        param.color = Leo::Color(0.3, 0.8, 0.1);

        m_func_ptr = std::make_unique<StdEvolution>(std::ref(m_core_ref), std::cref(param));
    }
};

}

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore2d<Ursa::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Ursa::create_window_properties("std lyapunov orbit"), argc, argv);

    window.show();
    window.run();
    return 0;
}
