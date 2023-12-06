///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <random>

#include "tools/test_tools.hpp"
#include <capd_utils/projection_map.hpp>
#include <capd_utils/timemap_wrapper.hpp>
#include <capd_utils/concat.hpp>
#include "pcr3bp_basic/levi_civita_coordinate_change.hpp"
#include "pcr3bp_basic/standard_system.hpp"
#include "pcr3bp_basic/regularized_system.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Compare flow value at point
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<size_t mu_index>
static void compare_with_standard_flow_at_point()
{
    using ScalarType = CapdUtils::Interval;
    using VectorType = CapdUtils::IVector;
    using MapT = CapdUtils::IMap;

    Pcr3bp::SetupParameters<MapT> setup(ScalarType(5) / ScalarType(100));
    MapT vf_reg = Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(mu_index, setup, false);
    MapT vf_std = Pcr3bp::StandardSystem<MapT>::createPositiveVectorField(setup);
    MapT gamma = LeviCivitaCoordinateChange<MapT>::create(mu_index, setup, true, false, false);
    MapT hamiltonian = Pcr3bp::StandardSystem<MapT>::createHamiltonian(setup);
    MapT R5toR4 = CapdUtils::ProjectionMap<MapT>::create(5, {0, 1, 2, 3});
    
    auto f = [&](VectorType reg_arg, ScalarType precision)
    {
        const VectorType std_arg = gamma(reg_arg);
        const ScalarType h = hamiltonian(std_arg)[0];
        const VectorType reg_arg5 = CapdUtils::Concat<MapT>::concat_vectors({ reg_arg, VectorType{ h } });

        const ScalarType u = reg_arg[0];
        const ScalarType v = reg_arg[1];
        const ScalarType factor = 4 * (u*u + v*v);

        const VectorType dU = vf_reg(reg_arg5);
        const VectorType diff = vf_std(std_arg) - gamma[reg_arg] * R5toR4(dU) / factor;

        CapdUtils::MaxNorm<MapT> norm;
        EXPECT_LT(norm(diff), precision);
    };

    std::mt19937_64 rng {};
    std::normal_distribution<double> dist {};

    for (int i = 0; i < 100; ++i)
    {
        const double u = dist(rng);
        const double v = dist(rng);
        const double pu = dist(rng);
        const double pv = dist(rng);

        f({ u, v, pu, pv }, 1.2e-12);
    }
}

}

TEST(Pcr3bp_basic, point_compare_regular_flow_at_1_with_standard_flow)
{
    Ursa::compare_with_standard_flow_at_point<1>();
}

TEST(Pcr3bp_basic, point_compare_regular_flow_at_2_with_standard_flow)
{
    Ursa::compare_with_standard_flow_at_point<2>();
}

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Compare pcr3bp flow in regularized coordinates with the one from standard coordinates
//!
//! @param mu_index index of mass at which the regularization takes place
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<size_t mu_index>
static void compare_with_standard_flow(Ursa::TimemapTestParams param)
{
    using ScalarType = CapdUtils::Interval;
    using VectorType = CapdUtils::IVector;
    using MapT = CapdUtils::IMap;

    Pcr3bp::SetupParameters<MapT> setup(ScalarType(5) / ScalarType(100));
    const ScalarType mu1 = setup.get_mu(1);
    const ScalarType mu2 = setup.get_mu(2);

    MapT R6toR4 = CapdUtils::ProjectionMap<MapT>::create(6, {0, 1, 2, 3});
    MapT reg_change6 = LeviCivitaCoordinateChange<MapT>::create(mu_index, setup, true, true, true);
    MapT hamiltonian = Pcr3bp::StandardSystem<MapT>::createHamiltonian(setup);
    MapT vf_std = Pcr3bp::StandardSystem<MapT>::createPositiveVectorField(setup);

    // create regularized flow map depending on mu_index at which
    // the regularization takes place
    MapT vf_reg = Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(mu_index, setup, true);

    CapdUtils::TimemapWrapper<MapT> std_flow(vf_std);
    CapdUtils::TimemapWrapper<MapT> reg_flow(vf_reg, param.time);

    VectorType U0 = param.initial;
    const VectorType X0 = R6toR4(reg_change6(U0));

    const VectorType H0 = hamiltonian(X0); // adjust H value
    U0[4] = H0[0];

    VectorType U = U0;
    VectorType X = X0;

    // evolve map in regular coordinates
    U = reg_flow(U);

    // extract final time and convert it into steps in standard flow evolution
    const Real final_time = CapdUtils::scalar_cast<Real>( U[5] );

    // evolve map in standard coordinates
    std_flow.set_time(final_time);
    X = std_flow(X);

    // compare points in both evolutions
    CapdUtils::MaxNorm<MapT> norm;
    EXPECT_LT(norm(X - R6toR4(reg_change6(U))), param.precision);
}

}

TEST(Pcr3bp_basic, compare_regular_flow_at_1_with_standard_flow_1)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.5, 0.1, 0.1, 1.1, 0.0, 0.0 };
    param.time = 0.1;
    param.precision = 1.5e-14;
    Ursa::compare_with_standard_flow<1>(param);
}

TEST(Pcr3bp_basic, compare_regular_flow_at_1_with_standard_flow_2)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.1, 2.1, 0.1, 1.1, 0.0, 0.0 };
    param.time = 0.5;
    param.precision = 1.9e-12;
    Ursa::compare_with_standard_flow<1>(param);
}

TEST(Pcr3bp_basic, compare_regular_flow_at_1_with_standard_flow_3)
{
    Ursa::TimemapTestParams param;
    param.initial = { -1.1, 0.5, 0.5, -0.5, 0.0, 0.0 };
    param.time = 0.5;
    param.precision = 1.8e-10;
    Ursa::compare_with_standard_flow<1>(param);
}

TEST(Pcr3bp_basic, compare_regular_flow_at_1_with_standard_flow_4)
{
    Ursa::TimemapTestParams param;
    param.initial = {  0.0, 0.5, 1.5, -1.5, 0.0, 0.0 };
    param.time = 0.5;
    param.precision = 2.9e-13;
    Ursa::compare_with_standard_flow<1>(param);
}

TEST(Pcr3bp_basic, compare_regular_flow_at_1_with_standard_flow_5)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.5, 0.0, 0.5, 1.0, 0.0, 0.0 };
    param.time = 0.5;
    param.precision = 7.0e-13;
    Ursa::compare_with_standard_flow<1>(param);
}

TEST(Pcr3bp_basic, compare_regular_flow_at_2_with_standard_flow_1)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.5, 0.1, 0.1, 1.1, 0.0, 0.0 };
    param.time = 0.1;
    param.precision = 3.2e-14;
    Ursa::compare_with_standard_flow<2>(param);
}

TEST(Pcr3bp_basic, compare_regular_flow_at_2_with_standard_flow_2)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.1, 2.1, 0.1, 1.1, 0.0, 0.0 };
    param.time = 0.5;
    param.precision = 4.4e-12;
    Ursa::compare_with_standard_flow<2>(param);
}

TEST(Pcr3bp_basic, compare_regular_flow_at_2_with_standard_flow_3)
{
    Ursa::TimemapTestParams param;
    param.initial = { -1.1, 0.5, 0.5, -0.5, 0.0, 0.0 };
    param.time = 0.5;
    param.precision = 6.9e-12;
    Ursa::compare_with_standard_flow<2>(param);
}

TEST(Pcr3bp_basic, compare_regular_flow_at_2_with_standard_flow_4)
{
    Ursa::TimemapTestParams param;
    param.initial = {  0.0, 0.5, 1.5, -1.5, 0.0, 0.0 };
    param.time = 0.5;
    param.precision = 4.6e-12;
    Ursa::compare_with_standard_flow<2>(param);
}

TEST(Pcr3bp_basic, compare_regular_flow_at_2_with_standard_flow_5)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.5, 0.0, 0.5, 1.0, 0.0, 0.0 };
    param.time = 0.5;
    param.precision = 8.5e-11;
    Ursa::compare_with_standard_flow<2>(param);
}

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Evolve pcr3bp in reg. coordinates forward time and then backward and compare final point with initial one
//!
//! @param mu_index index of mass at which the regularization takes place
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<size_t mu_index>
void evolve_forward_and_backward(Ursa::TimemapTestParams param)
{
    using ScalarType = CapdUtils::Interval;
    using VectorType = CapdUtils::IVector;
    using MapT = CapdUtils::IMap;

    Pcr3bp::SetupParameters<MapT> setup(ScalarType(5) / ScalarType(100));
    const ScalarType mu1 = setup.get_mu(1);
    const ScalarType mu2 = setup.get_mu(2);

    MapT R6toR4 = CapdUtils::ProjectionMap<MapT>::create(6, {0, 1, 2, 3});
    MapT reg_change6 = LeviCivitaCoordinateChange<MapT>::create(mu_index, setup, true, true, true);
    MapT hamiltonian = Pcr3bp::StandardSystem<MapT>::createHamiltonian(setup);

    MapT positive_vector_field = Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(mu_index, setup, false);
    MapT negative_vector_field = Pcr3bp::RegularizedSystem<MapT>::createNegativeVectorField(mu_index, setup, false);

    CapdUtils::TimemapWrapper<MapT> positive_flow(positive_vector_field, param.time);
    CapdUtils::TimemapWrapper<MapT> negative_flow(negative_vector_field, param.time);

    VectorType U0 = param.initial;
    const VectorType X0 = R6toR4(reg_change6(U0));

    const VectorType H0 = hamiltonian(X0); // adjust H value
    U0[4] = H0[0];

    const VectorType U1 = positive_flow(U0);
    const VectorType U2 = negative_flow(U1);

    CapdUtils::MaxNorm<MapT> norm;
    EXPECT_LT(norm(U2 - U0), param.precision);
}

}

TEST(Pcr3bp_basic, reg_flow_at_1_forward_and_backward_evolution)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.5, 0.0, 0.5, 1.0, 0.0 };
    param.time = 0.1;
    param.precision = 7.8e-15;
    Ursa::evolve_forward_and_backward<1>(param);
}

TEST(Pcr3bp_basic, reg_flow_at_2_forward_and_backward_evolution)
{
    Ursa::TimemapTestParams param;
    param.initial = { -0.5, 0.1, 0.1, 1.1, 0.0 };
    param.time = 0.1;
    param.precision = 1.3e-14;
    Ursa::evolve_forward_and_backward<2>(param);
}
