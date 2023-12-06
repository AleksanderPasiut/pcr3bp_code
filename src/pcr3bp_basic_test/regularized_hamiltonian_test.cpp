///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"
#include "pcr3bp_basic/levi_civita_coordinate_change.hpp"
#include "pcr3bp_basic/standard_system.hpp"
#include "pcr3bp_basic/regularized_system.hpp"

#include <capd_utils/projection_map.hpp>

namespace Ursa
{

template<typename MapT, size_t mu_index>
class Pcr3bpRegularizedHamiltonianTest
{
private:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;

    MapT m_hamiltonian;
    MapT m_R5toR4;
    MapT m_reg_change5;

    MapT m_ghv_condition_equation;

public:
    Pcr3bpRegularizedHamiltonianTest(const Pcr3bp::SetupParameters<MapT>& setup)
        : m_hamiltonian( Pcr3bp::StandardSystem<MapT>::createHamiltonian(setup) )
        , m_R5toR4( CapdUtils::ProjectionMap<MapT>::create(5, {0, 1, 2, 3}) )
        , m_reg_change5( LeviCivitaCoordinateChange<MapT>::create(mu_index, setup, true, true, false) )
        , m_ghv_condition_equation( Pcr3bp::RegularizedSystem<MapT>::createHamiltonian(mu_index, setup) )
    {}

    void check(const VectorType& U0, ScalarType precision)
    {
        const VectorType H0 = m_hamiltonian(m_R5toR4(m_reg_change5(U0)));
        const VectorType U1 = { U0[0], U0[1], U0[2], U0[3], H0[0] };

        CapdUtils::MaxNorm<MapT> norm;
        const VectorType V = m_ghv_condition_equation(U1);
        EXPECT_LT( norm( V ), precision );
    }
};

}

TEST(Pcr3bp_basic, reg_hamiltonian_test_1)
{
    using namespace Ursa;

    Pcr3bp::SetupParameters<IMap> setup;
    Pcr3bpRegularizedHamiltonianTest<IMap, 1> test(setup);

    test.check( IVector{ 0.1, 0.2, 0.3, 0.1 }, 1.4e-14 );
    test.check( IVector{ 0.5, -0.4, 0.7, 0.2 }, 8.0e-15 );
    test.check( IVector{ 0.3, -0.6, 0.1, 1.2 }, 1.1e-14 );
    test.check( IVector{ -0.8, -0.1, 1.1, 1.2 }, 1.1e-14 );
}

TEST(Pcr3bp_basic, reg_hamiltonian_test_2)
{
    using namespace Ursa;

    Pcr3bp::SetupParameters<IMap> setup(0.4);
    Pcr3bpRegularizedHamiltonianTest<IMap, 2> test(setup);

    test.check( IVector{ 0.1, 0.2, 0.3, 0.1 }, 8.9e-15 );
    test.check( IVector{ 0.5, -0.4, 0.7, 0.2 }, 5.8e-15 );
    test.check( IVector{ 0.3, -0.6, 0.1, 1.2 }, 6.3e-15 );
    test.check( IVector{ -0.8, -0.1, 1.1, 1.2 }, 1.3e-14 );
}
