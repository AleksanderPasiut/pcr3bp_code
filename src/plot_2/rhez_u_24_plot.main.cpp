///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_4d.hpp>

#include "plot_common/window_properties.hpp"

#include "rhez_u_24_plot.hpp"

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore4d<Pcr3bpProof::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Pcr3bpProof::create_window_properties("RHEZ U 24"),
        argc,
        argv,
        Leo::Color::White(),
        Serpent::Grid3d::Properties{
            .ruler_x = Leo::Ruler<>(-2, 2, 9, 1),
            .ruler_y = Leo::Ruler<>(-2, 2, 9, 1),
            .ruler_z = Leo::Ruler<>(-2, 2, 9, 1),
            .grid_color = Leo::Color::DarkGray(),
            .grid_width = 0.005f,
            .label_x = L".",
            .label_y = L".",
            .label_z = L".",
            .label_properties{
                .back_color = Leo::Color::White(),
                .fore_color = Leo::Color::DarkGray(),
                .halfheight = 0.05f
            }
        },
        Serpent::Colortray::get_properties_1d(
            Leo::Color::Profile::RedPurpleBlue11,
            Leo::Color::White(),
            Leo::Color::DarkGray(),
            Lyra::Colortray::get_default_placement_1d(),
            Leo::Ruler<double>(-1.0, +1.0, 3, 1))
    );

    window.show();
    window.run();
    return 0;
}
