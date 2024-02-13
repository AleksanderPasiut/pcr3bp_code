///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "covering_relations_test.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Part of interval arithmetic validation of the proof of Theorem 8, where we assert the existence of covering
//!        relations N_4 => N_5 => ... => N_K
//!
//!        In this computation we also perform a part of interval arithmetic validation of the proof of Lemma 8, where we
//!        assert that the trajectories shadowing the covering relations described above, do not intersect with the collision
//!        manifold.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(Pcr3bp_proof, homoclinic_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_homoclinic_coverings();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Part of interval arithmetic validation of the proof of Theorem 8, where we assert the existence of covering
//!        relations N_0 => N_1 => N_2
//!
//!        In this computation we also perform a part of interval arithmetic validation of the proof of Lemma 8, where we
//!        assert that the trajectories shadowing the covering relations described above, do not intersect with the collision
//!        manifold outside of h-set N_0.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(Pcr3bp_proof, periodic_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_periodic_coverings();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Part of interval arithmetic validation of the proof of Theorem 8, where we assert the existence of covering
//!        relation N_3 <= N_4
//!
//!        In this computation we also perform a part of interval arithmetic validation of the proof of Lemma 8, where we
//!        assert that the trajectories shadowing the covering relation described above, do not intersect with the collision
//!        manifold.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(Pcr3bp_proof, jump_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_jump_coverings();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Interval arithmetic validation of the proof of Lemma 11, where we assert the existence of covering relation
//!        SR_1 <= N_4
//!
//!        In this computation we also perform a part of interval arithmetic validation of the proof of Lemma 8, where we
//!        assert that the trajectories shadowing the covering relation described above, do not intersect with the collision
//!        manifold.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(Pcr3bp_proof, parallelogram_coverings_beginning)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.parallelogram_covering_beginning_check();
}
