///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "covering_relations_test.parallelogram_covering_derivative_check.hpp"

TEST(Pcr3bp_proof, parallelogram_coverings_derivative)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest_ParallelogramCoveringDerivativeCheck<IMap> test { setup };
    test.parallelogram_covering_derivative_check();
}
