///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_2d.hpp>

TEST(Plot, RegEvolutionRaw)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        Aquila::ParamConfig("u0", -5.0, 5.0, 0.0, 0.001, 9),
        Aquila::ParamConfig("v0", -5.0, 5.0, 0.0, 0.001, 9),
        Aquila::ParamConfig("pu0", -10.0, 10.0, 0.0, 0.001, 9),
        Aquila::ParamConfig("pv0", -10.0, 10.0, 0.3984, 0.001, 9),
        Aquila::ParamConfig("h", -10.0, 10.0, 0.0, 0.001, 9),
        Aquila::ParamConfig("t", 0.0, 100.0, 5.0, 0.1, 3),
        Aquila::ParamConfig("point_count", 0.0, 50000.0, 500.0, 1.0, 0)
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window2d window(paramset_config, "./plot_reg_evolution_raw");

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
