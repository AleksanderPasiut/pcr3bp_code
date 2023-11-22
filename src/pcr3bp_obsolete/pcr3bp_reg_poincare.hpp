///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/map_base.hpp>
#include <carina/poincare_wrapper.hpp>

#include <pcr3bp_basic/regularized_system.hpp>

#include <tools/direction.hpp>

namespace Ursa
{

template<typename MapT>
class Pcr3bpRegPoincareBase : Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpRegPoincareBase(const MapT& vf_reg, unsigned solver_order)
        : m_vf_reg(vf_reg)
        , m_poincare_v( m_vf_reg, solver_order, Carina::CoordinateSection<MapT>(m_vf_reg.dimension(), 1, ScalarType(0.0) ) )
    {}

    VectorType operator() (const VectorType& vec, MatrixType& mat)
    {
        this->assert_vector_size(vec, 5, "Pcr3bpRegPoincareBase vec vector size mismatch!");

        mat = MatrixType(5, 5);
        const VectorType ret = m_poincare_v(vec, mat);
        this->assert_matrix_size(mat, 5, 5, "Pcr3bpRegPoincareBase mat matrix size mismatch!");

        return ret;
    }

    VectorType operator() (const VectorType& vec)
    {
        this->assert_vector_size(vec, 5, "Pcr3bpRegPoincareBase vec vector size mismatch (2)!");

        const VectorType ret = m_poincare_v(vec);
        return ret;
    }

    ScalarType get_return_time(const VectorType& vec)
    {
        this->assert_vector_size(vec, 5, "Pcr3bpRegPoincareBase vec vector size mismatch (3)!");
        m_poincare_v(vec);
        return m_poincare_v.get_last_evaluation_return_time();
    }

    unsigned dimension() const noexcept { return 5; }
    unsigned imageDimension() const noexcept { return 5; }

private:
    MapT m_vf_reg;
    Carina::PoincareWrapper<MapT, Carina::CoordinateSection<MapT>> m_poincare_v;
};

template<typename MapT, Direction direction>
class Pcr3bpRegPoincare
{};

template<typename MapT>
class Pcr3bpRegPoincare<MapT, Direction::Positive> : public Pcr3bpRegPoincareBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpRegPoincare(const Pcr3bp::SetupParameters<MapT>& setup, unsigned solver_order = 20)
        : Pcr3bpRegPoincareBase<MapT>(
            Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(2, setup, false),
            solver_order)
    {}
};

template<typename MapT>
class Pcr3bpRegPoincare<MapT, Direction::Negative> : public Pcr3bpRegPoincareBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpRegPoincare(const Pcr3bp::SetupParameters<MapT>& setup, unsigned solver_order = 20)
        : Pcr3bpRegPoincareBase<MapT>(
            Pcr3bp::RegularizedSystem<MapT>::createNegativeVectorField(2, setup, false),
            solver_order)
    {}
};

template <typename MapT>
using Pcr3bpRegPoincarePositiveU = Pcr3bpRegPoincare<MapT, Direction::Positive>;

template <typename MapT>
using Pcr3bpRegPoincareNegativeU = Pcr3bpRegPoincare<MapT, Direction::Negative>;

}
