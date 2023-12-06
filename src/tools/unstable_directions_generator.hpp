///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/gauss.hpp>

#include "power_iteration.hpp"

namespace Carina
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Generate unstable positive and unstable negative directions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, typename MapPosT, typename MapNegT>
class UnstableDirectionsGenerator
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    UnstableDirectionsGenerator(
        MapPosT map_pos,
        MapNegT map_neg,
        VectorType initial_dir,
        unsigned power_iteration_steps = 50)
    {
        assert_with_exception(map_pos.dimension() == map_neg.dimension());
        assert_with_exception(map_pos.dimension() == map_pos.imageDimension());
        assert_with_exception(map_neg.dimension() == map_neg.imageDimension());

        MatrixType m_pos(map_pos.imageDimension(), map_pos.dimension());
        const VectorType p_pos = map_pos( VectorType(map_pos.dimension()), m_pos );

        MatrixType m_neg(map_neg.imageDimension(), map_neg.dimension());
        const VectorType p_neg = map_neg( VectorType(map_neg.dimension()), m_neg );
        
        const MatrixType m_pos_total = gaussInverseMatrix<MapT>(m_neg) * m_pos;
        const MatrixType m_neg_total = gaussInverseMatrix<MapT>(m_pos) * m_neg;

        m_unstable_dir_pos = PowerIteration<MapT>::evaluate( m_pos_total, initial_dir, power_iteration_steps);
        m_expansion_pos_factor = (m_pos_total * m_unstable_dir_pos).euclNorm();

        m_unstable_dir_neg = PowerIteration<MapT>::evaluate( m_neg_total, initial_dir, power_iteration_steps);
        m_expansion_neg_factor = (m_neg_total * m_unstable_dir_neg).euclNorm();
    }

    VectorType get_unstable_pos() const
    {
        return m_unstable_dir_pos;
    }

    VectorType get_unstable_neg() const
    {
        return -m_unstable_dir_neg;
    }

    ScalarType get_expansion_pos_factor() const
    {
        return m_expansion_pos_factor;
    }

    ScalarType get_expansion_neg_factor() const
    {
        return m_expansion_neg_factor;
    }

private:
    VectorType m_unstable_dir_pos;
    VectorType m_unstable_dir_neg;

    ScalarType m_expansion_pos_factor;
    ScalarType m_expansion_neg_factor;

};

}
