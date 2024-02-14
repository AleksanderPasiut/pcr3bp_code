///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Check parallelogram covering conditions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class ParallelogramCoveringChecker
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    ParallelogramCoveringChecker(MatrixType aligned_poincare_map_derivative)
    {
        const MatrixType& der = aligned_poincare_map_derivative;

        const ScalarType alpha = 5.09;
        const ScalarType beta = 0.195;
        const ScalarType rho = 0.197;
        const ScalarType c = 0.021;

        EXPECT_TRUE( 0 < beta );
        EXPECT_TRUE( beta < rho );
        EXPECT_TRUE( alpha > 2*c + rho );
        EXPECT_TRUE( c + rho < 1 );

        EXPECT_TRUE( der(1,1) > alpha );
        EXPECT_TRUE( der(1,2) < 0 );
        EXPECT_TRUE( der(1,2) > -c );
        EXPECT_TRUE( der(2,1) > 0 );
        EXPECT_TRUE( der(2,1) < c );
        EXPECT_TRUE( der(2,2) < rho );
        EXPECT_TRUE( der(2,2) > beta );
    }
};

}

