///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_2d.hpp>

TEST(Plot_Test, Gen2dPlot)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        Aquila::ParamConfig("plot_type", 0.0, 20.0, 7.0, 1.0, 0),
        Aquila::ParamConfig("point_count", 0.0, 50000.0, 100.0, 1.0, 0),

        Aquila::ParamConfig("thickness", 0.0, 1.0, 0.002, 0.0001, 9),

        Aquila::ParamConfig("span_x1 (pow10)", -100.0, 100.0, -13.0, 1.0, 10),
        Aquila::ParamConfig("span_x2 (pow10)", -100.0, 100.0, -10.5, 1.0, 10),
        Aquila::ParamConfig("scale (pow10)", -100.0, 100.0, -10.0, 1.0, 10),

        Aquila::ParamConfig("shift_x1", -10.0, 10.0, 0.0, 1.0, 3),
        Aquila::ParamConfig("shift_x2", -10.0, 10.0, 0.0, 1.0, 3),

        Aquila::ParamConfig("s min", -3.0, 3.0, -5.0e-9, 1.0e-9, 9),
        Aquila::ParamConfig("s max", -3.0, 3.0, +5.0e-9, 1.0e-9, 9),
        Aquila::ParamConfig("t min", -3.0, 3.0, -1.0e-5, 1.0e-9, 9),
        Aquila::ParamConfig("t max", -3.0, 3.0, +1.0e-5, 1.0e-9, 9),

        Aquila::ParamConfig("alpha1", -3.2, 3.2,  1.40156, 1.0e-9, 9),
        Aquila::ParamConfig("alpha2", -3.2, 3.2, -1.40156, 1.0e-9, 9)
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window2d window(paramset_config, "./gen2dplot");

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
