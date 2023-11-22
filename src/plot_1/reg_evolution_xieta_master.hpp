///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_2d.hpp>

TEST(Plot, RegEvolutionXieta)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        Aquila::ParamConfig("xi", -2.5, 2.5, 1.0, 0.001, 9),
        Aquila::ParamConfig("eta", -2.5, 2.5, 0.0, 0.001, 9),
        Aquila::ParamConfig("point_count", 0.0, 50000.0, 200.0, 1.0, 0),
        Aquila::ParamConfig("dual poincare", 0.0, 1.0, 0.0, 1.0, 0),
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window2d window(paramset_config, "./plot_reg_evolution_xieta");

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
