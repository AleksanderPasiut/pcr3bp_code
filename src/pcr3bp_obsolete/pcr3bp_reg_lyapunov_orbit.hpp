///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/map_base.hpp>
#include <capd_utils/timemap_wrapper.hpp>
#include <capd_utils/poincare_wrapper.hpp>
#include <capd_utils/enp_map.hpp>
#include <capd_utils/pne_map.hpp>
#include <capd_utils/identity_map.hpp>

#include <pcr3bp_basic/regularized_system.hpp>

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (u0, v0, pu0, pv0, u1, v1, pu1, pv1, u2, v2, pu2, pv2, h0) -> (timemap(U0) - U1, poincare(U1) - U2)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LyapunovOrbitRegBase : CapdUtils::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LyapunovOrbitRegBase(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType t)
        : m_vf_reg( Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(2, setup, false) )
        , m_timemap_u( m_vf_reg, t, 60)
        , m_poincare_v( m_vf_reg, 60, CapdUtils::CoordinateSection<MapT>(m_vf_reg.dimension(), 1, ScalarType(0.0)) )
        , m_timemap_u_pne(13,
            CapdUtils::IdxList<size_t>{ 0, 1, 2, 3, 12 },
            CapdUtils::IdxList<int>::create(10, 0, 5),
            std::ref(m_timemap_u))
        , m_poincare_v_pne(13,
            CapdUtils::IdxList<size_t>{ 4, 5, 6, 7, 12 },
            CapdUtils::IdxList<int>::create(10, 5, 5),
            std::ref(m_poincare_v))
        , m_id_1_pne(13,
            CapdUtils::IdxList<size_t>{ 4, 5, 6, 7, 12 },
            CapdUtils::IdxList<int>::create(10, 0, 5),
            5)
        , m_id_2_pne(13,
            CapdUtils::IdxList<size_t>{ 8, 9, 10, 11, 12 },
            CapdUtils::IdxList<int>::create(10, 5, 5),
            5)
    {}

    VectorType operator() (const VectorType& vec, MatrixType& mat)
    {
        this->assert_vector_size(vec, 13, "LyapunovOrbitRegBase vec vector size mismatch!");

        MatrixType m1(10, 13);
        const VectorType v1 = m_timemap_u_pne(vec, m1);
        this->assert_matrix_size(m1, 10, 13, "LyapunovOrbitRegBase m1 matrix size mismatch!");

        MatrixType m2(10, 13);
        const VectorType v2 = m_poincare_v_pne(vec, m2);
        this->assert_matrix_size(m2, 10, 13, "LyapunovOrbitRegBase m2 matrix size mismatch!");

        MatrixType m3(10, 13);
        const VectorType v3 = m_id_1_pne(vec, m3);
        this->assert_matrix_size(m3, 10, 13, "LyapunovOrbitRegBase m3 matrix size mismatch!");

        MatrixType m4(10, 13);
        const VectorType v4 = m_id_2_pne(vec, m4);
        this->assert_matrix_size(m4, 10, 13, "LyapunovOrbitRegBase m4 matrix size mismatch!");

        mat = m1 + m2 - m3 - m4;
        return v1 + v2 - v3 - v4;
    }

    VectorType operator() (const VectorType& vec)
    {
        MatrixType dummy;
        return (*this)(vec, dummy);
    }

    unsigned dimension() const noexcept
    {
        return 13;
    }

    unsigned imageDimension() const noexcept
    {
        return 10;
    }

private:
    MapT m_vf_reg;

    CapdUtils::TimemapWrapper<MapT> m_timemap_u;
    CapdUtils::PoincareWrapper<MapT, CapdUtils::CoordinateSection<MapT>> m_poincare_v;

    CapdUtils::PNE<MapT, decltype(m_timemap_u)&> m_timemap_u_pne;
    CapdUtils::PNE<MapT, decltype(m_poincare_v)&> m_poincare_v_pne;

    CapdUtils::PNE<MapT, CapdUtils::IdentityMap<MapT>> m_id_1_pne;
    CapdUtils::PNE<MapT, CapdUtils::IdentityMap<MapT>> m_id_2_pne;
};

template<typename MapT>
class LyapunovOrbitRegCollision : public CapdUtils::ENP<MapT, LyapunovOrbitRegBase<MapT>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LyapunovOrbitRegCollision(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType pv0, ScalarType t = 1.0)
        : CapdUtils::ENP<MapT, LyapunovOrbitRegBase<MapT>>(
            VectorType{ zero(), zero(), zero(), pv0, zero(), zero(), zero(), zero(), zero(), zero(), zero(), zero(), zero() },
            { -1, -1, -1, -1,
              0,  1,  2,  3,
              -1, -1, -1, -1, 4 },
            { 0, 1, 2, 3, 7 },

            std::cref(setup), t)
    {}

private:
    static constexpr ScalarType zero()
    {
        return ScalarType(0.0);
    }
};

}
