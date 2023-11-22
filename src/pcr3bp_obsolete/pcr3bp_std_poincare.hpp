///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/map_base.hpp>
#include <carina/poincare_wrapper.hpp>

#include <pcr3bp_basic/standard_system.hpp>

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Std. pcr3bp poincare positive map on y == 0 section
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Pcr3bpStdPoincare : Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bpStdPoincare(const Pcr3bp::SetupParameters<MapT>& setup)
        : m_vf_std( Pcr3bp::StandardSystem<MapT>::createPositiveVectorField(setup) )
        , m_poincare( m_vf_std, 40, Carina::CoordinateSection<MapT>(m_vf_std.dimension(), 1, ScalarType(0.0) ) )
    {}

    VectorType operator() (const VectorType& vec)
    {
        this->assert_vector_size(vec, 4, "Pcr3bpStdPoincare vec vector size mismatch (1)!");

        const VectorType v1 = m_poincare(vec);
        return v1;
    }

    VectorType operator() (const VectorType& vec, MatrixType& mat)
    {
        this->assert_vector_size(vec, 4, "Pcr3bpStdPoincare vec vector size mismatch (2)!");

        const VectorType v1 = m_poincare(vec, mat);
        this->assert_matrix_size(mat, 4, 4, "Pcr3bpStdPoincare m1 matrix size mismatch!");

        return v1;
    }

    ScalarType get_return_time(const VectorType& vec)
    {
        this->assert_vector_size(vec, 4, "Pcr3bpStdPoincare vec vector size mismatch (3)!");
        m_poincare(vec);
        return m_poincare.get_last_evaluation_return_time();
    }

    unsigned dimension() const noexcept
    {
        return m_vf_std.dimension();
    }

    unsigned imageDimension() const noexcept
    {
        return m_vf_std.imageDimension();
    }

private:
    MapT m_vf_std;
    Carina::PoincareWrapper<MapT, Carina::CoordinateSection<MapT>> m_poincare;
};

}
