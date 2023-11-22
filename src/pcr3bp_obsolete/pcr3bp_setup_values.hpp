///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"
#include <carina/readable_interval.hpp>

#include <pcr3bp_basic/setup_parameters.hpp>
#include <unordered_map>


namespace std
{

template<>
struct hash<Ursa::Interval>
{
    size_t operator() (const Ursa::Interval& arg) const
    {
        std::hash<double> ha {};        
        return ha(arg.leftBound()) + ha(arg.rightBound());
    }
};

}

namespace Ursa
{

template<typename ScalarType>
inline ScalarType compute_sqrt(const ScalarType& arg);

template<>
inline Real compute_sqrt(const Real& arg)
{
    return std::sqrt(arg);
}

template<>
inline Interval compute_sqrt(const Interval& interval)
{
    return capd::intervals::sqrt(interval);
}

template<typename MapT>
class Pcr3bpSetupValues : public Pcr3bp::SetupParameters<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;

    Pcr3bpSetupValues(ScalarType mu1 = Pcr3bp::SetupParameters<MapT>::get_default_mu1() )
        : Pcr3bp::SetupParameters<MapT>(mu1)
    {
        {
            const ScalarType mu = ScalarType(1.0) / ScalarType(100);
            m_collision_orbit_h0[mu] = Carina::scalar_cast<ScalarType, Interval>(Carina::ReadableInterval<Interval>("bfe6b201d51bd5b6", "bfe6b201d51bd164"));

            IVector positive_dir(2);
            positive_dir[0] = Carina::ReadableInterval<Interval>("bfc58f3772f48905", "bfc58f12680f8890");
            positive_dir[1] = Carina::ReadableInterval<Interval>("bfef8af7bf28d4af", "bfef8af62a11aa11");
            m_positive_unstable_dir[mu] = Carina::vector_cast<VectorType>(positive_dir);

            IVector negative_dir(2);
            negative_dir[0] = Carina::ReadableInterval<Interval>("bfc58f377cf6a561", "bfc58f125e0c3b6c");
            negative_dir[1] = Carina::ReadableInterval<Interval>("3fef8af629a436e9", "3fef8af7bf965415");
            m_negative_unstable_dir[mu] = Carina::vector_cast<VectorType>(negative_dir);
        }

        {
            const ScalarType mu = 1.0 / 100;
            m_collision_orbit_h0[mu] = Carina::scalar_cast<ScalarType, Interval>(Carina::ReadableInterval<Interval>("bfe6b201d51bd5b6", "bfe6b201d51bd164"));

            IVector positive_dir(2);
            positive_dir[0] = Carina::ReadableInterval<Interval>("bfc58f3772f48905", "bfc58f12680f8890");
            positive_dir[1] = Carina::ReadableInterval<Interval>("bfef8af7bf28d4af", "bfef8af62a11aa11");
            m_positive_unstable_dir[mu] = Carina::vector_cast<VectorType>(positive_dir);

            IVector negative_dir(2);
            negative_dir[0] = Carina::ReadableInterval<Interval>("bfc58f377cf6a561", "bfc58f125e0c3b6c");
            negative_dir[1] = Carina::ReadableInterval<Interval>("3fef8af629a436e9", "3fef8af7bf965415");
            m_negative_unstable_dir[mu] = Carina::vector_cast<VectorType>(negative_dir);
        }

        {
            const ScalarType mu = ScalarType(1.0) / ScalarType(82);
            m_collision_orbit_h0[mu] = Carina::scalar_cast<ScalarType, Interval>(Carina::ReadableInterval<Interval>("bfe6c0fe27bfd27d", "bfe6c0fe27bfcddf"));
            m_initial_collision_orbit_h0[mu] = Carina::scalar_cast<ScalarType, Interval>(Carina::ReadableInterval<Interval>("bfe6c0fe27bfd27b", "bfe6c0fe27bfcde0"));

            IVector positive_dir(2);
            positive_dir[0] = Carina::ReadableInterval<Interval>("bfc546639b8b4bc9", "bfc5458b5131a954");
            positive_dir[1] = Carina::ReadableInterval<Interval>("bfef8e164f8595e3", "bfef8e0d3277e76f");
            m_positive_unstable_dir[mu] = Carina::vector_cast<VectorType>(positive_dir);

            IVector negative_dir(2);
            negative_dir[0] = Carina::ReadableInterval<Interval>("bfc546640b0d04a9", "bfc5458ae1b08d04");
            negative_dir[1] = Carina::ReadableInterval<Interval>("3fef8e0d2dc50589", "3fef8e1654383ff5");
            m_negative_unstable_dir[mu] = Carina::vector_cast<VectorType>(negative_dir);
        }

        {
            const ScalarType mu = 1.0 / 82;
            m_collision_orbit_h0[mu] = Carina::scalar_cast<ScalarType, Interval>(Carina::ReadableInterval<Interval>("bfe6c0fe27bfd27d", "bfe6c0fe27bfcddf"));
            m_initial_collision_orbit_h0[mu] = Carina::scalar_cast<ScalarType, Interval>(Carina::ReadableInterval<Interval>("bfe6c0fe27bfd27a", "bfe6c0fe27bfcddd"));

            IVector positive_dir(2);
            positive_dir[0] = Carina::ReadableInterval<Interval>("bfc54663a4e010a9", "bfc5458b47dd9300");
            positive_dir[1] = Carina::ReadableInterval<Interval>("bfef8e164fea3357", "bfef8e0d32133e84");
            m_positive_unstable_dir[mu] = Carina::vector_cast<VectorType>(positive_dir);

            IVector negative_dir(2);
            negative_dir[0] = Carina::ReadableInterval<Interval>("bfc546642b3c3085", "bfc5458ac1817fe0");
            negative_dir[1] = Carina::ReadableInterval<Interval>("3fef8e0d2c69d532", "3fef8e16559360c2");
            m_negative_unstable_dir[mu] = Carina::vector_cast<VectorType>(negative_dir);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Get energy value of the Lyapunov orbit passing through collision
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ScalarType get_collision_orbit_h0() const
    {
        return m_collision_orbit_h0.at( this->get_mu(1) );
    }

    ScalarType get_initial_collision_orbit_h0() const
    {
        return m_initial_collision_orbit_h0.at( this->get_mu(1) );
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Get value of regularized momentum pv component of the Lyapunov orbit passing through collision
    //!
    //! @details This value is applicable for regularization performed at mu2 mass.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ScalarType get_collision_pv() const
    {
        return compute_sqrt(this->get_mu(2) * 8);
    }

    VectorType get_positive_unstable_dir() const
    {
        return m_positive_unstable_dir.at( this->get_mu(1) );
    }

    VectorType get_negative_unstable_dir() const
    {
        return m_negative_unstable_dir.at( this->get_mu(1) );
    }

private:
    std::unordered_map<ScalarType, ScalarType> m_collision_orbit_h0;
    std::unordered_map<ScalarType, ScalarType> m_initial_collision_orbit_h0;

    std::unordered_map<ScalarType, VectorType> m_positive_unstable_dir;
    std::unordered_map<ScalarType, VectorType> m_negative_unstable_dir;
};

}
