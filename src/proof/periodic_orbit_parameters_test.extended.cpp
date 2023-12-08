///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "periodic_orbit_parameters_test.hpp"

TEST(Pcr3bp_intermediate, periodic_orbit_parameters_test_nonrigorous)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();
    LyapunovOrbitRegCollisionSetup<RMap> setup( 58696.0 / 65536 );
}

TEST(Extended_Pcr3bp_intermediate, periodic_orbit_parameters_test_rigorous_2)
{
    using namespace Pcr3bpProof;

    LyapunovOrbitRegCollisionSetup<IMap> setup( 7.0 / 8, false );
}

TEST(Exteded_Pcr3bp_intermediate, periodic_orbit_parameters_test_rigorous_3)
{
    using namespace Pcr3bpProof;

    LyapunovOrbitRegCollisionSetup<IMap> setup( 1.0 / 2, false );
}
