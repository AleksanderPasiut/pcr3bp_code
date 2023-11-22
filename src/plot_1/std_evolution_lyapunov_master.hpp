///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_2d.hpp>

TEST(Plot, StdEvolutionLyapunov)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        Aquila::ParamConfig("x0", -1.0, 1.0, 0.9, 0.001, 9),
        Aquila::ParamConfig("point_count", 0.0, 50000.0, 20.0, 1.0, 0),
        Aquila::ParamConfig("steps", 0.0, 1000.0, 5.0, 1.0, 0)
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window2d window(paramset_config, "./plot_std_evolution_lyapunov");

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
