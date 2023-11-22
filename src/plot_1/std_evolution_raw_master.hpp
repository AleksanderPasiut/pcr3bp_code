///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_2d.hpp>

TEST(Plot, StdEvolutionRaw)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        Aquila::ParamConfig("x0", -2.0, 2.0, 0.5, 0.001, 9),
        Aquila::ParamConfig("y0", -2.0, 2.0, 0.0, 0.001, 9),
        Aquila::ParamConfig("px0", -100.0, 100.0, 0.0, 0.001, 9),
        Aquila::ParamConfig("py0", -100.0, 100.0, 1.0, 0.001, 9),
        Aquila::ParamConfig("t", 0.0, 1000.0, 1.0, 0.1, 3),
        Aquila::ParamConfig("point_count", 0.0, 50000.0, 10.0, 1.0, 0)
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window2d window(paramset_config, "./plot_std_evolution_raw");

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
