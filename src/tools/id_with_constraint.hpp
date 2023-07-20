///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/projection_map.hpp>
#include <carina/constrained_function.hpp>

namespace Carina
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @todo
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, typename ConstraintT>
class IdWithConstraint : public MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    IdWithConstraint(ConstraintT constraint) : m_constraint(constraint)
    {}

    VectorType operator() (const VectorType& vec) override
    {
        this->assert_vector_size(vec, m_internal.dimension(), "IdWithConstraint vec vector size mismatch! (1)");
        return m_internal(vec);
    }

    VectorType operator() (const VectorType& vec, MatrixType& mat) override
    {
        this->assert_vector_size(vec, m_internal.dimension(), "IdWithConstraint vec vector size mismatch! (2)");
        return m_internal(vec, mat);
    }

    unsigned dimension() const override
    {
        return m_internal.dimension();
    }

    unsigned imageDimension() const override
    {
        return m_internal.imageDimension();
    }

private:
    ConstraintT m_constraint;

    MapT m_id { ProjectionMap<MapT>::create(m_constraint.dimension(), IdxList<size_t>::create(0, m_constraint.dimension())) };

    ConstrainedFunction<MapT, MapT&, ConstraintT> m_internal
    {
        std::ref(m_id),
        std::ref(m_constraint),
        VectorType(m_constraint.imageDimension())
    };
};

}
