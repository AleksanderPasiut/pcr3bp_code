///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_2d.hpp>

TEST(Plot, RegEvolutionLyapunovMulti)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        Aquila::ParamConfig("u0", -1.0, 1.0, 0.0, 0.001, 9),
        Aquila::ParamConfig("point_count", 0.0, 50000.0, 800.0, 1.0, 0),
        Aquila::ParamConfig("steps", 0.0, 1000.0, 5.0, 1.0, 0),
        Aquila::ParamConfig("option", 0.0, 1000.0, 0.0, 1.0, 0),
        Aquila::ParamConfig("point size", 0.0, 1.0, 0.02, 0.001, 3),
        Aquila::ParamConfig("evolution time", 0.0, 10.0, 6.272, 0.001, 3),
        Aquila::ParamConfig("selected point", -1.0, 50.0, -1.0, 1.0, 0)
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window2d window(paramset_config, "./plot_reg_evolution_lyapunov_multi");

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
