///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"
#include <carina/capd/norm.hpp>
#include <carina/timemap_wrapper.hpp>
#include "pcr3bp_basic/standard_system.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Check if hamiltonian value remains constant along evolution
//! calculated with pcr3bp positive flow map with standard coordinates
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
static void check_hamiltonian_constantness(Ursa::TimemapTestParams param)
{
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;

    const Pcr3bp::SetupParameters<MapT> setup(0.05);

    MapT hamiltonian = Pcr3bp::StandardSystem<MapT>::createHamiltonian(setup);

    MapT vector_field = Pcr3bp::StandardSystem<MapT>::createPositiveVectorField(setup);

    Carina::TimemapWrapper<MapT> flow(vector_field, param.time);

    Carina::MaxNorm<MapT> norm;

    const VectorType X0 = param.initial;
    const VectorType H0 = hamiltonian(X0);

    const VectorType X1 = flow(X0);
    const VectorType H1 = hamiltonian(X1);

    EXPECT_LT(norm(H1 - H0), param.precision);

    const VectorType X2 = flow(X0);
    const VectorType H2 = hamiltonian(X2);

    EXPECT_LT(norm(H2 - H0), param.precision);
}

}

TEST(Pcr3bp_basic, std_flow_hamiltonian_constantness_1)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.5, 0.1, 0.1, 1.1 };
    param.time = 0.1;
    param.precision = 1.3e-14;
    Ursa::check_hamiltonian_constantness<Carina::IMap>(param);
}
TEST(Pcr3bp_basic, std_flow_hamiltonian_constantness_2)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.1, 1.1, 0.5, 0.3 };
    param.time = 0.5;
    param.precision = 1.5e-14;
    Ursa::check_hamiltonian_constantness<Carina::IMap>(param);
}
TEST(Pcr3bp_basic, std_flow_hamiltonian_constantness_3)
{
    Ursa::TimemapTestParams param;
    param.initial = { -1.1, 0.5, -0.5, 0.1 };
    param.time = 0.2;
    param.precision = 5.0e-15;
    Ursa::check_hamiltonian_constantness<Carina::IMap>(param);
}
TEST(Pcr3bp_basic, std_flow_hamiltonian_constantness_4)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.1, -0.8, 1.5, 1.3 };
    param.time = 0.1;
    param.precision = 1.9e-14;
    Ursa::check_hamiltonian_constantness<Carina::IMap>(param);
}

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Evolve pcr3bp forward in time and then backward and compare
//! final point with initial one
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
static void check_forward_and_backward_evolution(Ursa::TimemapTestParams param)
{
    const Pcr3bp::SetupParameters<MapT> setup(0.05);

    IMap positive_vector_field = Pcr3bp::StandardSystem<MapT>::createPositiveVectorField(setup);
    IMap negative_vector_field = Pcr3bp::StandardSystem<MapT>::createNegativeVectorField(setup);

    Carina::TimemapWrapper<IMap> positive_flow(positive_vector_field, param.time);
    Carina::TimemapWrapper<IMap> negative_flow(negative_vector_field, param.time);

    const IVector X0 = param.initial;
    const IVector X1 = positive_flow(X0);
    const IVector X2 = negative_flow(X1);

    Carina::MaxNorm<IMap> norm;
    EXPECT_LT(norm(X2 - X0), param.precision);
}

}

TEST(Pcr3bp_basic, std_flow_forward_and_backward_evolution_1)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.5, 0.1, 0.1, 1.1 };
    param.time = 0.1;
    param.precision = 5.9e-15;
    Ursa::check_forward_and_backward_evolution<Ursa::IMap>(param);
}
TEST(Pcr3bp_basic, std_flow_forward_and_backward_evolution_2)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.1, 1.1, 0.5, 0.3 };
    param.time = 0.5;
    param.precision = 1.4e-14;
    Ursa::check_forward_and_backward_evolution<Ursa::IMap>(param);
}
TEST(Pcr3bp_basic, std_flow_forward_and_backward_evolution_3)
{
    Ursa::TimemapTestParams param;
    param.initial = { -1.1, 0.5, -0.5, 0.1 };
    param.time = 0.2;
    param.precision = 3.8e-15;
    Ursa::check_forward_and_backward_evolution<Ursa::IMap>(param);
}
TEST(Pcr3bp_basic, std_flow_forward_and_backward_evolution_4)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.1, -0.8, 1.5, 1.3 };
    param.time = 0.1;
    param.precision = 6.5e-15;
    Ursa::check_forward_and_backward_evolution<Ursa::IMap>(param);
}
