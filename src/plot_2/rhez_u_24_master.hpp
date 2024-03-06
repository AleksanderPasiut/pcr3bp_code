///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_4d.hpp>

TEST(Plot, Rhez_u_24)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        Aquila::ParamConfig("point thickness", 0.0, 0.1, 0.05, 0.0001, 4),
        Aquila::ParamConfig("reg evolution point count", 0.0, 100000.0, 2000.0, 1.0, 0),
        Aquila::ParamConfig("reg evolution thickness", 0.0, 0.05, 0.005, 0.0001, 4),


        Aquila::ParamConfig("select short path section", -1.0, 3.0, -1.0, 1.0, 0),
        Aquila::ParamConfig("select long path section", -1.0, 27.0, -1.0, 1.0, 0),

        Aquila::ParamConfig("select short path section CE", -1.0, 3.0, -1.0, 1.0, 0),
        Aquila::ParamConfig("select long path section CE", -1.0, 27.0, -1.0, 1.0, 0),

        Aquila::ParamConfig("section span", 0.0001, 1000.0, 1.0, 0.0001, 4),
        // Aquila::ParamConfig("section scale", 0.0001, 1000.0, 1.0, 0.0001, 4),

        Aquila::ParamConfig("evolution time", 0.0, 10.0, 2.636, 0.0001, 6),

        Aquila::ParamConfig("scale", -20.0, 20.0, 0.0, 0.01, 2),

        Aquila::ParamConfig("show periodic orbit evo.", 0.0, 1.0, 0.0, 1.0, 0),
        Aquila::ParamConfig("show homoclinic orbit evo.", 0.0, 1.0, 0.0, 1.0, 0),

        Aquila::ParamConfig("show periodic orbit origins", 0.0, 1.0, 0.0, 1.0, 0),
        Aquila::ParamConfig("show homoclinic orbit origins", 0.0, 1.0, 0.0, 1.0, 0),

        Aquila::ParamConfig("highlight periodic orbit point", -1.0, 100.0, -1.0, 1.0, 0),
        Aquila::ParamConfig("highlight homoclinic orbit point", -1.0, 100.0, -1.0, 1.0, 0),

        Aquila::ParamConfig("centerpoint index", -1.0, 100.0, -1.0, 1.0, 0)
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window4d window(paramset_config, "./plot_rhez_u_24", 75);

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
