///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/composite_map.hpp>
#include <capd_utils/affine_map.hpp>
#include <capd_utils/image_sum.hpp>
#include <capd_utils/local_poincare_wrapper.hpp>

#include "id_with_constraint.hpp"

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @todo
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LocalPoincare4_Constraint : public CapdUtils::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LocalPoincare4_Constraint(
        MapT& constraint,
        const CapdUtils::LocalCoordinateSystem<MapT>& src_coordsys)
            : m_constraint(constraint)
            , m_src_coordsys(src_coordsys)
    {
        assert_with_exception(m_constraint.dimension() == 4);
        assert_with_exception(m_constraint.imageDimension() == 1);
        assert_with_exception(m_src_coordsys.get_origin().dimension() == 4);
        assert_with_exception(m_src_coordsys.get_origin().dimension() == 4);
    }

    VectorType operator() (const VectorType& vec) override
    {
        return m_extension_to_4(vec);
    }

    VectorType operator() (const VectorType& vec, MatrixType& der) override
    {
        return m_extension_to_4(vec, der);
    }

    unsigned dimension() const override
    {
        return m_extension_to_4.dimension();
    }

    unsigned imageDimension() const override
    {
        return m_extension_to_4.imageDimension();
    }

private:
    MapT& m_constraint;
    
    const CapdUtils::LocalCoordinateSystem<MapT> m_src_coordsys;

    CapdUtils::AffineMap<MapT> m_linear_src { m_src_coordsys };

    CapdUtils::AffineMap2<MapT> m_affine_section_constraint
    {
        m_src_coordsys.get_origin(),
        CapdUtils::Extract<MapT>::get_vvector( m_src_coordsys.get_directions_matrix(), 3 )
    };

    CapdUtils::CompositeMap<MapT,
        decltype(m_linear_src)&,
        decltype(m_affine_section_constraint)&> m_affine_section_constraint_src
    {
        std::ref(m_linear_src),
        std::ref(m_affine_section_constraint)
    };

    CapdUtils::CompositeMap<MapT, decltype(m_linear_src)&, MapT&> m_constraint_src
    {
        std::ref(m_linear_src),
        std::ref(m_constraint)
    };

    CapdUtils::ImageSum<MapT, decltype(m_affine_section_constraint_src)&, decltype(m_constraint_src)&> m_dual_constraint
    {
        std::ref(m_affine_section_constraint_src),
        std::ref(m_constraint_src)
    };

    CapdUtils::IdWithConstraint<MapT, decltype(m_dual_constraint)&> m_extension_to_4
    {
        std::ref(m_dual_constraint)
    };
};

}
