///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <gtest/gtest.h>
#include <taurus/window_4d.hpp>

TEST(Plot, Rhez_u_34)
{
    Aquila::ParamsetConfig<double> paramset_config
    {
        // Aquila::ParamConfig("u_span", 0.0, 1.0, 0.997, 0.001, 9)
    };

    Glib::RefPtr<Gtk::Application> app =
        Gtk::Application::create("org.ampasiut.ursa_plot");

    Taurus::Window4d window(paramset_config, "./plot_rhez_u_34", 50);

    const int ret = app->run(window);
    EXPECT_EQ(0, ret);
}
