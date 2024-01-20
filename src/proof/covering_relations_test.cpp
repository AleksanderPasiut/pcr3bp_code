///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "covering_relations_test.hpp"

TEST(Pcr3bp_proof, homoclinic_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_homoclinic_coverings();
}

TEST(Pcr3bp_proof, periodic_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_periodic_coverings();
}

TEST(Pcr3bp_proof, jump_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_jump_coverings();
}

TEST(Pcr3bp_proof, parallelogram_coverings_derivative)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.parallelogram_covering_derivative_check();
}

TEST(Pcr3bp_proof, parallelogram_coverings_beginning)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.parallelogram_covering_beginning_check();
}

TEST(Pcr3bp_proof, parallelogram_coverings_collision)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.collision_manifold_derivative_check();
}
