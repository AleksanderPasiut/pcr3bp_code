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
#include "objects/coordinate_systems_origins.hpp"
#include "objects/homoclinic_orbit.hpp"

#include "proof/periodic_orbit_parameters.hpp"

namespace Ursa
{

class CoreInterior : CoreInteriorBase
{
public:
    CoreInterior(Lyra::Core2d& core_ref)
        : CoreInteriorBase()
        , m_core_ref(core_ref)
        , m_masses(core_ref, m_setup, 0.03f, Leo::Color(0.1, 0.1, 0.4), Leo::Color(0.6, 0.0, 0.0))
    {}

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        CoreInteriorBase::set_param(packet_vector);

        m_evolutions.clear();

        const size_t point_count = static_cast<size_t>(this->get_param(1));
        const size_t steps = static_cast<size_t>(this->get_param(2));
        const size_t option = static_cast<size_t>(this->get_param(3));
        const double point_size = this->get_param(4);
        const double evolution_time = this->get_param(5);

        std::vector<double> u0_vec = { -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15 };

        for (const double& u0 : u0_vec)
        {
            if (option == 0)
            {
                RegEvolutionParam param = get_reg_evolution_param(u0, steps, point_count);
                m_evolutions.emplace_back(std::ref(m_core_ref), std::cref(param));
            }
            else
            {
                m_evolutions.clear();
            }

            if (option == 1)
            {
                RegEvolutionParam param = get_reg_evolution_param(u0, steps, point_count);
                m_evolutions_std.emplace_back(std::ref(m_core_ref), std::cref(param));
            }
            else
            {
                m_evolutions_std.clear();
            }

            if (option == 2)
            {
                m_coordinate_systems_origins = std::make_unique<CoordinateSystemsOrigins>( std::ref(m_core_ref), m_setup, point_size );
            }
            else
            {
                m_coordinate_systems_origins.reset();
            }
        }

        if (option == 3)
        {
            m_homoclinic_orbit = std::make_unique<HomoclinicOrbit>(std::ref(m_core_ref), evolution_time, point_count);
        }
        else
        {
            m_homoclinic_orbit.reset();
        }
    }

private:
    RegEvolutionParam get_reg_evolution_param(double u0, size_t steps, size_t point_count)
    {
        const RVector PV = LyapunovOrbitRegParam::calculate(m_setup, u0, steps);

        RegEvolutionParam param;
        param.setup = m_setup;
        param.u0 = PV[0];
        param.v0 = PV[1];
        param.pu0 = PV[2];
        param.pv0 = PV[3];
        param.h = PV[4];

        Pcr3bpRegPoincarePositiveU<RMap> poincare( m_setup );
        Real t = poincare.get_return_time(RVector{ param.u0, param.v0, param.pu0, param.pv0, param.h });
        param.t = 2 * t;

        param.point_count = point_count;

        param.point_thickness = 0.0f;//5e-3f;
        param.line_thickness = 0.006f;
        param.point_subcount = 10;
        param.color = (u0 == 0.0) ? Leo::Color(1.0, 0.0, 0.0) : Leo::Color(0.3, 0.1, 0.8);
        return param;
    }

    Lyra::Core2d& m_core_ref;

    Pcr3bp::SetupParameters<RMap> m_setup {};
    RegMasses m_masses;

    std::list<RegEvolution> m_evolutions {};
    std::list<RegEvolutionWithCoordChange> m_evolutions_std {};

    std::unique_ptr<CoordinateSystemsOrigins> m_coordinate_systems_origins {};

    std::unique_ptr<HomoclinicOrbit> m_homoclinic_orbit {};

};

}

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore2d<Ursa::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Ursa::create_window_properties("reg lyapunov orbit multi"),
        argc,
        argv,
        Leo::Color(1.0, 1.0, 1.0),
        Serpent::Grid2d::Properties{
            .ruler_x = Leo::Ruler<>(-1.5, 1.5, 7, 1),
            .ruler_y = Leo::Ruler<>(-1.5, 1.5, 7, 1),
            .grid_color = Leo::Color(0.1, 0.1, 0.1),
            .grid_width = 0.005f,
            .label_x = L"u",
            .label_y = L"v",
            .label_properties{
                .back_color = Leo::Color(1.0, 1.0, 1.0),
                .fore_color = Leo::Color(0.1, 0.1, 0.1),
                .halfheight = 0.03f,
                .align = Serpent::Align::Center
            }
        });

    window.show();
    window.run();
    return 0;
}
