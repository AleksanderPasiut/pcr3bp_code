///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "local_poincare4_projection_base.hpp"
#include "auxiliary_functions.hpp"

#include <capd_utils/local_coordinate_system.hpp>
#include <capd_utils/affine_map.hpp>
#include <capd_utils/composite_map.hpp>

namespace Pcr3bpProof
{

template<typename MapT>
class LocalPoincare4_Projection_Spec : public LocalPoincare4_ProjectionBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LocalPoincare4_Projection_Spec(const CapdUtils::LocalCoordinateSystem<MapT>& dst_coordsys) : m_dst_coordsys_4_dim(dst_coordsys)
    {}

    virtual VectorType operator() (const VectorType& vec) override
    {
        print_var(vec);
        print_var( m_affine_map(vec) );

        return m_internal(vec);
    }

    virtual VectorType operator() (const VectorType& vec, MatrixType& der) override
    {
        print_var(vec);
        print_var( m_affine_map(vec) );

        return m_internal(vec, der);
    }

    virtual unsigned dimension() const override
    {
        return m_internal.dimension();
    }

    virtual unsigned imageDimension() const override
    {
        return m_internal.imageDimension();
    }

private:
    CapdUtils::LocalCoordinateSystem<MapT> m_dst_coordsys_4_dim;

    CapdUtils::AffineMap<MapT> m_affine_map
    {
        m_dst_coordsys_4_dim.get_origin(),
        m_dst_coordsys_4_dim.get_directions_matrix()
    };

    MapT m_constraint_inverse
    {
        AuxiliaryFunctions<MapT>::create_psi0_inverse( Psi0_Coefficients<MapT>::get().get_d_coeffs() )
    };

    CapdUtils::CompositeMap<MapT,
        decltype(m_affine_map)&,
        decltype(m_constraint_inverse)&> m_internal
    {
        std::ref(m_affine_map),
        std::ref(m_constraint_inverse)
    };
};

}
