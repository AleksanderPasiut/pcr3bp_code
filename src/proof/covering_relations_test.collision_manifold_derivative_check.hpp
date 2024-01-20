///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"

#include "tools/auxiliary_functions.hpp"
#include "tools/local_poincare4_constraint_spec.hpp"

#include "covering_relations_test_base.hpp"
#include "covering_relation_checker.hpp"

namespace Pcr3bpProof
{

template<typename MapT>
class CoveringRelationsTest : public CoveringRelationsTestBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = CapdUtils::LocalCoordinateSystem<MapT>;

    CoveringRelationsTest(const CoveringRelationsSetup& setup)
        : CoveringRelationsTestBase<MapT>(setup)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check collision manifold derivative
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void collision_manifold_derivative_check()
    {
        const ScalarType L = this->m_basic_objects.m_parallelogram_coverings_parameters.L;

        MapT eta = AuxiliaryFunctions<MapT>::eta( L );

        LocalPoincare4_Constraint_Spec<MapT> map_E0
        {
            std::ref(this->m_basic_objects.m_hamiltonian_reg2),
            std::ref( this->m_periodic_orbit_coordsys.at(0) )
        };

        MatrixType out(4, 2);
        map_E0( VectorType(2), out );

        print_var(out);

        CapdUtils::AffineMap<MapT> map_L0 { this->m_periodic_orbit_coordsys.at(0) };

        MapT& map_C = this->m_basic_objects.m_collision_condition;

        CapdUtils::CompositeMap<MapT,
            decltype(eta)&,
            decltype(map_E0)&,
            decltype(map_L0)&,
            decltype(map_C)&> map_chi(
                std::ref(eta),
                std::ref(map_E0),
                std::ref(map_L0),
                std::ref(map_C)
            );

        MatrixType der(3, 2);
        map_chi(N * this->m_gain_factor, der);
        const ScalarType fx = der(3,1);
        const ScalarType fy = der(3,2);

        EXPECT_TRUE( fx < 0 );
        EXPECT_TRUE( fy < 0 );

        print_var ( fx );
        print_var ( fy );
    }
};

}
