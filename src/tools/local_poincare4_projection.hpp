///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "local_poincare4_projection_base.hpp"

#include <capd_utils/projection_map.hpp>

namespace Pcr3bpProof
{

template<typename MapT>
class LocalPoincare4_Projection : public LocalPoincare4_ProjectionBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LocalPoincare4_Projection(const CapdUtils::LocalCoordinateSystem<MapT>& dst_coordsys)
    {}

    virtual VectorType operator() (const VectorType& vec) override
    {
        return m_projection_to_2(vec);
    }

    virtual VectorType operator() (const VectorType& vec, MatrixType& der) override
    {
        return m_projection_to_2(vec, der);
    }

    virtual unsigned dimension() const override
    {
        return m_projection_to_2.dimension();
    }

    virtual unsigned imageDimension() const override
    {
        return m_projection_to_2.imageDimension();
    }

private:
    MapT m_projection_to_2
    {
        CapdUtils::ProjectionMap<MapT>::create( 4, { 0, 1 } )
    };
};

}
