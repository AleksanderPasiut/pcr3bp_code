///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/enp_map.hpp>

#include <pcr3bp_obsolete/pcr3bp_rhez_common.hpp>
#include <pcr3bp_obsolete/pcr3bp_reg_poincare.hpp>

#include "pcr3bp_obsolete/pcr3bp_setup_values.hpp"
#include "pcr3bp_rhez_alpha.hpp"
#include "pcr3bp_rhez_xieta.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Reg. poincare map in specified coordinates
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, Direction direction, CoordType coord_type>
class Pcr3bpRegPoincareCoord : Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpRegPoincareCoord(const Pcr3bp::SetupParameters<MapT>& setup, unsigned solver_order = 20)
        : m_poincare(setup, solver_order)
        , m_rhez(setup)
        , m_rhez_inv(setup)
    {}

    VectorType operator() (const VectorType& vec, MatrixType& mat)
    {
        this->assert_vector_size(vec, 3, "Pcr3bpRegPoincareCoord vec vector size mismatch!");

        MatrixType m1(5, 3);
        const VectorType v1 = m_rhez(vec, m1);
        this->assert_matrix_size(m1, 5, 3, "Pcr3bpRegPoincareCoord m1 matrix size mismatch!");

        MatrixType m2(5, 5);
        const VectorType v2 = m_poincare(v1, m2);
        this->assert_matrix_size(m2, 5, 5, "Pcr3bpRegPoincareCoord m2 matrix size mismatch!");

        MatrixType m3(3, 5);
        const VectorType v3 = m_rhez_inv(v2, m3);
        this->assert_matrix_size(m3, 3, 5, "Pcr3bpRegPoincareCoord m3 matrix size mismatch!");

        mat = m3 * m2 * m1;

        return v3;
    }

    VectorType operator() (const VectorType& vec)
    {
        this->assert_vector_size(vec, 3, "Pcr3bpRegPoincareCoord vec vector size mismatch (2)!");

        const VectorType v1 = m_rhez(vec);
        const VectorType v2 = m_poincare(v1);
        const VectorType v3 = m_rhez_inv(v2);

        return v3;
    }

    unsigned dimension() const noexcept { return 3; }
    unsigned imageDimension() const noexcept { return 3; }

private:
    Pcr3bpRegPoincare<MapT, direction> m_poincare;

    Pcr3bpRhez<MapT, 2, coord_type> m_rhez;
    Pcr3bpRhezInverse<MapT, 2, coord_type> m_rhez_inv;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Reg. poincare map in specified coordinates for collision energy
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, Direction direction, CoordType coord_type>
class Pcr3bpRegPoincareCoord_CE : public Carina::ENP<MapT, Pcr3bpRegPoincareCoord<MapT, direction, coord_type>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpRegPoincareCoord_CE(const Pcr3bpSetupValues<MapT>& setup, unsigned solver_order = 20)
        : Carina::ENP<MapT, Pcr3bpRegPoincareCoord<MapT, direction, coord_type>>(
            VectorType{ ScalarType(0.0), ScalarType(0.0), setup.get_collision_orbit_h0() },
            { 0, 1, -1 },
            { 0, 1 },

        std::cref(setup), solver_order)
    {}
};

template<typename MapT>
using Pcr3bpRegPoincarePositiveU_Alpha = Pcr3bpRegPoincareCoord<MapT, Direction::Positive, CoordType::UAlpha>;

template<typename MapT>
using Pcr3bpRegPoincareNegativeU_Alpha = Pcr3bpRegPoincareCoord<MapT, Direction::Negative, CoordType::UAlpha>;

template<typename MapT>
using Pcr3bpRegPoincarePositiveU_Alpha_CE = Pcr3bpRegPoincareCoord_CE<MapT, Direction::Positive, CoordType::UAlpha>;

template<typename MapT>
using Pcr3bpRegPoincareNegativeU_Alpha_CE = Pcr3bpRegPoincareCoord_CE<MapT, Direction::Negative, CoordType::UAlpha>;

template<typename MapT>
using Pcr3bpRegPoincarePositiveU_XiEta = Pcr3bpRegPoincareCoord<MapT, Direction::Positive, CoordType::UXiEta>;

template<typename MapT>
using Pcr3bpRegPoincareNegativeU_XiEta = Pcr3bpRegPoincareCoord<MapT, Direction::Negative, CoordType::UXiEta>;

template<typename MapT>
using Pcr3bpRegPoincarePositiveU_XiEta_CE = Pcr3bpRegPoincareCoord_CE<MapT, Direction::Positive, CoordType::UXiEta>;

template<typename MapT>
using Pcr3bpRegPoincareNegativeU_XiEta_CE = Pcr3bpRegPoincareCoord_CE<MapT, Direction::Negative, CoordType::UXiEta>;

}
