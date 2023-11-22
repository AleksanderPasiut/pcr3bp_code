///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_3d.hpp>

TEST(Plot, Rhez_u)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        Aquila::ParamConfig("h", -2.0, 2.0, -0.709229389, 0.01, 9),
        Aquila::ParamConfig("u_span", 0.0, 1.0, 0.999, 0.001, 9),
        Aquila::ParamConfig("manifold_select", 0.0, 65535.0, 2.0, 1.0, 0),

        Aquila::ParamConfig("level", 0.0, 100.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("points", 0.0, 10000.0, 10.0, 1.0, 0),
        Aquila::ParamConfig("sub points", 0.0, 10000.0, 3.0, 1.0, 0),
        Aquila::ParamConfig("min multiplier", -30.0, 30.0, -30.0, 0.1, 3),
        Aquila::ParamConfig("max multiplier", -30.0, 30.0, -6.0, 0.1, 3),

        Aquila::ParamConfig("show interp. nodes", 0.0, 1.0, 1.0, 1.0, 0),

        Aquila::ParamConfig("x", -3.0, 3.0, 0.0, 0.000001, 9),
        Aquila::ParamConfig("y", -3.0, 3.0, 2.814249, 0.00001, 9),
        Aquila::ParamConfig("z", -3.0, 3.0, 0.0, 0.00001, 9),

        Aquila::ParamConfig("range min", -10.0, +10.0, -0.1, 0.01, 2),
        Aquila::ParamConfig("range max", -10.0, +10.0, +1.1, 0.01, 2),

        Aquila::ParamConfig("s min", -3.0, 3.0, 0.0, 0.000001, 9),
        Aquila::ParamConfig("s max", -3.0, 3.0, 0.0, 0.000001, 9),
        Aquila::ParamConfig("t min", -3.0, 3.0, 0.0, 0.000001, 9),
        Aquila::ParamConfig("t max", -3.0, 3.0, 0.0, 0.000001, 9)
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window3d window(paramset_config, "./plot_rhez_u");

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
