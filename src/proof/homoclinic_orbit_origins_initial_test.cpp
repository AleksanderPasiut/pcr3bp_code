///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"
#include "tools/print_bootstrap.hpp"
#include "tools/variable_printer.hpp"
#include "homoclinic_orbit_origins_initial_generator.hpp"
#include "homoclinic_orbit_origins_initial.hpp"

TEST(Pcr3bp_intermediate, homoclinic_orbit_origins_initial_test)
{
    using namespace Pcr3bpProof;

    using MapT = RMap;
    using ScalarType = MapT::ScalarType;

    capd::rounding::DoubleRounding::roundNearest();

    PeriodicOrbitCoordsysGenerator<MapT> setup {};
    
    HomoclinicOrbitOriginsInitialGenerator<MapT> homoclinic_orbit_origins_initial_generator { setup.get_coordsys_container() };

    #if 1

    std::ofstream ostr("homoclinic_orbit_parameters.cpp.generated");

    if (ostr)
    {
        const auto& homoclinic_orbit_points = homoclinic_orbit_origins_initial_generator.get_points();
        for (auto it = homoclinic_orbit_points.begin(); ;)
        {
            CapdUtils::BootstrapPrint<MapT>::print(ostr, *it);

            ++it;
            if (it != homoclinic_orbit_points.end())
            {
                ostr << ",\n";
            }
            else
            {
                ostr << "\n";
                break;
            }
        }

        ostr.close();
    }

    CapdUtils::VariablePrinter<MapT>::print(
            "homoclinic_orbit_total_expansion_factor.txt",
            "Total expansion factor along homoclinic orbit",
            homoclinic_orbit_origins_initial_generator.get_total_expansion_factor());

    for (auto v : homoclinic_orbit_origins_initial_generator.get_points())
    {
        print_var(v);
    }

    #endif

    HomoclinicOrbitOriginsInitial<MapT> homoclinic_orbit_origins_initial {};

    EXPECT_EQ(
        homoclinic_orbit_origins_initial.get_total_expansion_factor(),
        homoclinic_orbit_origins_initial_generator.get_total_expansion_factor() );

    ASSERT_EQ( 
        homoclinic_orbit_origins_initial.get_points().size(),
        homoclinic_orbit_origins_initial_generator.get_points().size() );

    for (size_t i = 0; i < homoclinic_orbit_origins_initial.get_points().size(); ++i)
    {
        EXPECT_EQ(
            homoclinic_orbit_origins_initial.get_points().at(i),
            homoclinic_orbit_origins_initial_generator.get_points().at(i)) << i;
    }
}
