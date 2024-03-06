///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_2d.hpp>

TEST(Plot, Coverings)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        Aquila::ParamConfig("thickness", 0.0, 1.0, 0.002, 0.0001, 9),
        Aquila::ParamConfig("multiplier", -100.0, 100.0, 0.0, 1.0, 3),

        Aquila::ParamConfig("center x", -10.0, 10.0, 0.0, 0.001, 9),
        Aquila::ParamConfig("center y", -10.0, 10.0, 0.0, 0.001, 9),
        Aquila::ParamConfig("magnifier", -100.0, 100.0, 0.0, 1.0, 3),

        Aquila::ParamConfig("center_id", -1.0, 100.0, -1.0, 1.0, 3),

        Aquila::ParamConfig("show L0", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L1", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L2", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L3", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L4", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L5", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L6", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L7", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L8", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L9", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L10", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L11", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L12", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L13", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L14", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L15", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L16", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L17", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show L18", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show P0", 0.0, 1.0, 1.0, 1.0, 0),
        Aquila::ParamConfig("show P1", 0.0, 1.0, 1.0, 1.0, 0)
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window2d window(paramset_config, "./plot_coverings");

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
