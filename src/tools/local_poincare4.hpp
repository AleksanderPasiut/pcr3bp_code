///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"

#include <capd_utils/poincare_wrapper.hpp>
#include <capd_utils/timemap_wrapper.hpp>
#include <capd_utils/extension_map.hpp>
#include <capd_utils/projection_map.hpp>
#include <capd_utils/constrained_function.hpp>

#include "id_with_constraint.hpp"
#include "affine_poincare_map.hpp"

#include "local_poincare4_constraint.hpp"
#include "local_poincare4_constraint_spec.hpp"

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Implementation of local Poincare evaluation 
//! @details Implementation of the function:
//!
//!     \psi_m^{-1} \circ P \circ \psi_k
//!
//! where indices k and m are implicitly specified with source and destination coordinate systems.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LocalPoincare4 : public CapdUtils::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LocalPoincare4(
        MapT& vector_field,
        MapT& constraint,
        unsigned order,
        const CapdUtils::LocalCoordinateSystem<MapT>& src_coordsys,
        const CapdUtils::LocalCoordinateSystem<MapT>& dst_coordsys,
        bool specialized)
            : m_vector_field(vector_field)
            , m_constraint(constraint)
            , m_order(order)
            , m_src_coordsys(src_coordsys)
            , m_dst_coordsys(dst_coordsys)
            , m_specialized(specialized)
    {
        assert_with_exception(m_vector_field.dimension() == 4);
        assert_with_exception(m_vector_field.imageDimension() == 4);
        assert_with_exception(m_constraint.dimension() == 4);
        assert_with_exception(m_constraint.imageDimension() == 1);
        assert_with_exception(m_src_coordsys.get_origin().dimension() == 4);
        assert_with_exception(m_src_coordsys.get_origin().dimension() == 4);
    }

    VectorType operator() (const VectorType& vec) override
    {
        return m_affine_poincare_2(vec);
    }

    VectorType operator() (const VectorType& vec, MatrixType& der) override
    {
        return m_affine_poincare_2(vec, der);
    }

    unsigned dimension() const override
    {
        return m_extension_to_4.dimension();
    }

    unsigned imageDimension() const override
    {
        return m_projection_to_2.imageDimension();
    }

    ScalarType get_last_evaluation_return_time() const
    {
        return m_affine_poincare.get_last_evaluation_return_time();
    }

    void operator() (const VectorType& vec, ScalarType time, CapdUtils::SolutionCurve<MapT>& solution_curve)
    {
        const VectorType e = m_extension_to_4(vec);
        const VectorType v = m_affine_src(e);

        m_timemap.set_time(time);
        m_timemap(v, solution_curve);
    }

private:
    MapT& m_vector_field;
    MapT& m_constraint;

    unsigned m_order;
    
    const CapdUtils::LocalCoordinateSystem<MapT> m_src_coordsys;
    const CapdUtils::LocalCoordinateSystem<MapT> m_dst_coordsys;

    const bool m_specialized;

    CapdUtils::AffinePoincareMap<MapT> m_affine_poincare
    {
        m_vector_field,
        m_order,
        m_src_coordsys,
        m_dst_coordsys,
    };

    using LocalPoincare4_Constraint_BaseType = LocalPoincare4_Constraint_Base<MapT>;
    using LocalPoincare4_Constraint_BaseTypePtr = std::unique_ptr<LocalPoincare4_Constraint_BaseType>;
    using LocalPoincare4_Constraint_Type = LocalPoincare4_Constraint<MapT>;
    using LocalPoincare4_Constraint_SpecType = LocalPoincare4_Constraint_Spec<MapT>;

    LocalPoincare4_Constraint_BaseTypePtr m_extension_to_4_ptr
    {
        [this]() -> LocalPoincare4_Constraint_BaseTypePtr
        {
            return m_specialized ?
                LocalPoincare4_Constraint_BaseTypePtr(std::make_unique<LocalPoincare4_Constraint_SpecType>(
                    std::ref(m_constraint),
                    std::ref(m_src_coordsys)
                )) :
                LocalPoincare4_Constraint_BaseTypePtr(std::make_unique<LocalPoincare4_Constraint_Type>(
                    std::ref(m_constraint),
                    std::ref(m_src_coordsys)
                ));
        }()
    };

    LocalPoincare4_Constraint_BaseType& m_extension_to_4
    {
        *m_extension_to_4_ptr
    };

    MapT m_projection_to_2
    {
        CapdUtils::ProjectionMap<MapT>::create( 4, { 0, 1 } )
    };

    CapdUtils::CompositeMap<MapT,
        decltype(m_extension_to_4)&,
        decltype(m_affine_poincare)&,
        decltype(m_projection_to_2)&> m_affine_poincare_2
    {
        std::ref(m_extension_to_4),
        std::ref(m_affine_poincare),
        std::ref(m_projection_to_2)
    };

    CapdUtils::AffineMap<MapT> m_affine_src
    {
        m_src_coordsys
    };

    CapdUtils::TimemapWrapper<MapT> m_timemap
    {
        m_vector_field,
        0.0,
        m_order
    };
};

}
