///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <random>

#include "tools/test_tools.hpp"
#include "pcr3bp_basic/levi_civita_coordinate_change.hpp"
#include "pcr3bp_basic/levi_civita_inverse_coordinate_change.hpp"

namespace Ursa
{

class LeviCivitaCoordinateChangeTest
{
public:
    using MapT = IMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Perform coordinate change from standard to regularized and back and compare the results.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void check_coordinate_change_to_reg_and_back(const VectorType& X, ScalarType precision)
    {
        const VectorType U = m_rci(X);
        const VectorType X2 = m_rcf(U);
        EXPECT_LT(m_norm(X - X2), precision);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Perform coordinate change from regularized to standard and back and compare the results.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void check_coordinate_change_to_std_and_back(const VectorType& U, ScalarType precision)
    {
        const VectorType X = m_rcf(U);
        const VectorType U2 = m_rci(X);
        EXPECT_LT(m_norm(U - U2), precision);
    }

private:
    Pcr3bp::SetupParameters<MapT> m_setup { ScalarType(1) / ScalarType(10) };
    
    LeviCivitaInverseCoordinateChange<MapT> m_rci{ m_setup.get_x(2) };
    MapT m_rcf { LeviCivitaCoordinateChange<MapT>::create(2, m_setup, true, true, false) };

    CapdUtils::MaxNorm<MapT> m_norm {};
};

}

TEST(Pcr3bp_basic, coordinate_change_std_to_reg_and_back)
{
    using namespace Ursa;

    std::mt19937_64 rng {};
    std::normal_distribution<double> dist {};

    LeviCivitaCoordinateChangeTest test {};

    for (int i = 0; i < 100;)
    {
        const double x = dist(rng);
        const double y = dist(rng);

        if ((x+0.1)*(x+0.1) + y*y > 0.0025)
        {
            const double px = dist(rng);
            const double py = dist(rng);
            test.check_coordinate_change_to_reg_and_back({ x, y, px, py, 0.0 }, 4.1e-12);
            ++i;
        }
    }
}

TEST(Pcr3bp_basic, coordinate_change_reg_to_std_and_back)
{
    using namespace Ursa;

    std::mt19937_64 rng {};
    std::normal_distribution<double> dist {};

    LeviCivitaCoordinateChangeTest test {};

    for (int i = 0; i < 100;)
    {
        const double u = std::abs(dist(rng));
        const double v = dist(rng);

        if (u > 0.05)
        {
            const double pu = dist(rng);
            const double pv = dist(rng);
            test.check_coordinate_change_to_std_and_back({ u, v, pu, pv, 0.0 }, 3.2e-12);
            ++i;
        }
    }
}

TEST(Pcr3bp_basic, coordinate_change_test)
{
    using namespace Ursa;

    std::mt19937_64 rng {};
    std::normal_distribution<double> dist {};

    using MapT = IMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;

    Pcr3bp::SetupParameters<MapT> m_setup { ScalarType(1) / ScalarType(10) };

    MapT m_lcf { LeviCivitaCoordinateChange<MapT>::create(2, m_setup, true, false, false) };

    CapdUtils::MaxNorm<MapT> m_norm {};

    auto test = [m_norm, m_lcf](MapT m_lci, VectorType U, ScalarType precision)
    {
        const VectorType X = m_lcf(U);
        const VectorType U2 = m_lci(X);
        EXPECT_LT(m_norm(U - U2), precision);
    };

    {
        SCOPED_TRACE("PositiveU test");
        MapT m_lci { LeviCivitaCoordinateChange<MapT>::createInverse(2, m_setup, true, false, LeviCivitaCoordinateChangeInverseVariant::PositiveU ) };

        for (int i = 0; i < 100;)
        {
            const double u = std::abs( dist(rng) );
            const double v = dist(rng);

            if (u > 0.05)
            {
                const double pu = dist(rng);
                const double pv = dist(rng);
                const VectorType U = { u, v, pu, pv };

                test(m_lci, U, 7.5e-14);
                ++i;
            }
        }
    }

    {
        SCOPED_TRACE("NegativeU test");
        MapT m_lci { LeviCivitaCoordinateChange<MapT>::createInverse(2, m_setup, true, false, LeviCivitaCoordinateChangeInverseVariant::NegativeU ) };

        for (int i = 0; i < 100;)
        {
            const double u = -std::abs( dist(rng) );
            const double v = dist(rng);

            if (-u > 0.05)
            {
                const double pu = dist(rng);
                const double pv = dist(rng);
                const VectorType U = { u, v, pu, pv };

                test(m_lci, U, 7.5e-14);
                ++i;
            }
        }
    }

    {
        SCOPED_TRACE("PositiveV test");
        MapT m_lci { LeviCivitaCoordinateChange<MapT>::createInverse(2, m_setup, true, false, LeviCivitaCoordinateChangeInverseVariant::PositiveV ) };

        for (int i = 0; i < 100;)
        {
            const double v = std::abs( dist(rng) );
            const double u = dist(rng);

            if (v > 0.05)
            {
                const double pu = dist(rng);
                const double pv = dist(rng);
                const VectorType U = { u, v, pu, pv };

                test(m_lci, U, 2.9e-14);
                ++i;
            }
        }
    }

    {
        SCOPED_TRACE("NegativeV test");
        MapT m_lci { LeviCivitaCoordinateChange<MapT>::createInverse(2, m_setup, true, false, LeviCivitaCoordinateChangeInverseVariant::NegativeV ) };

        for (int i = 0; i < 100;)
        {
            const double v = -std::abs( dist(rng) );
            const double u = dist(rng);

            if (-v > 0.05)
            {
                const double pu = dist(rng);
                const double pv = dist(rng);
                const VectorType U = { u, v, pu, pv };

                test(m_lci, U, 2.7e-13);
                ++i;
            }
        }
    }
}
