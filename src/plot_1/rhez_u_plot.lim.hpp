///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/map_base.hpp>
#include <carina/newton_method.obsolete.hpp>
#include <carina/projection_map.hpp>

#include "tools/test_tools.hpp"
#include "capd_renderable.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief (s, t) -> (pu, pv, u)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Lim : Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Lim()
    {
        m_origin = VectorType{ 0.0, -3.1041492906675194341, 0.99796369890049541063 };

        const VectorType a = { 0.1933917702993516452,-3.4906330386670481047, 0.99831688272000340678 };
        const VectorType b = { -0.19339177029903714677,-3.4906330386665462839, 0.99831688272009277973 };

        m_matrix = MatrixType{
            { a[0] - m_origin[0], b[0] - m_origin[0] },
            { a[1] - m_origin[1], b[1] - m_origin[1] },
            { a[2] - m_origin[2], b[2] - m_origin[2] }
        };

        const ScalarType inv_det = m_matrix(1,1)*m_matrix(2,2) - m_matrix(1,2)*m_matrix(2,1);

        m_matrix_inv = MatrixType{
            { m_matrix(2,2), -m_matrix(1,2) },
            { -m_matrix(2,1), m_matrix(1,1) }
        };

        m_matrix_inv /= inv_det;
    }

    VectorType operator() (const VectorType& vec) const
    {
        this->assert_vector_size(vec, 2, "Lim vec vector size mismatch!");
        return m_matrix * vec + m_origin;
    }

    VectorType operator() (const VectorType& vec, MatrixType& mat) const
    {
        mat = m_matrix;
        return (*this)(vec);
    }

    VectorType compute_inverse(const VectorType& vec, MatrixType& mat) const
    {
        mat = m_matrix_inv;
        return m_matrix_inv * (vec - VectorType{ m_origin[0], m_origin[1] } );
    }

private:
    VectorType m_origin;
    MatrixType m_matrix;

    MatrixType m_matrix_inv;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief (pu, pv) -> (s, t)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LimInverse : public Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LimInverse(Lim<MapT>& lim) : m_lim(lim)
    {}

    VectorType operator() (const VectorType& vec)
    {
        this->assert_vector_size(vec, 2, "LimInverse vec vector size mismatch (1)!");

        MatrixType dummy;
        return m_lim.compute_inverse(vec, dummy);
    }

    VectorType operator() (const VectorType& vec, MatrixType& mat)
    {
        this->assert_vector_size(vec, 2, "LimInverse vec vector size mismatch (2)!");

        return m_lim.compute_inverse(vec, mat);
    }

private:
    Lim<MapT>& m_lim;
};

}
