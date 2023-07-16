///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/composite_map.hpp>
#include <carina/affine_map.hpp>
#include <carina/image_sum.hpp>
#include <carina/local_poincare_wrapper.hpp>

#include "id_with_constraint.hpp"

namespace Ursa
{

template<typename MapT>
class LocalPoincare4_Constraint : public Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LocalPoincare4_Constraint(
        MapT& constraint,
        const Carina::LocalCoordinateSystem<MapT>& src_coordsys)
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
    
    const Carina::LocalCoordinateSystem<MapT> m_src_coordsys;

    Carina::AffineMap<MapT> m_linear_src { m_src_coordsys };

    Carina::AffineMap2<MapT> m_affine_section_constraint
    {
        m_src_coordsys.get_origin(),
        Carina::Extract<MapT>::get_vvector( m_src_coordsys.get_directions_matrix(), 3 )
    };

    Carina::CompositeMap<MapT,
        decltype(m_linear_src)&,
        decltype(m_affine_section_constraint)&> m_affine_section_constraint_src
    {
        std::ref(m_linear_src),
        std::ref(m_affine_section_constraint)
    };

    Carina::CompositeMap<MapT, decltype(m_linear_src)&, MapT&> m_constraint_src
    {
        std::ref(m_linear_src),
        std::ref(m_constraint)
    };

    Carina::ImageSum<MapT, decltype(m_affine_section_constraint_src)&, decltype(m_constraint_src)&> m_dual_constraint
    {
        std::ref(m_affine_section_constraint_src),
        std::ref(m_constraint_src)
    };

    Carina::IdWithConstraint<MapT, decltype(m_dual_constraint)&> m_extension_to_4
    {
        std::ref(m_dual_constraint)
    };
};

}
