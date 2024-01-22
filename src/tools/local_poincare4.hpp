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

#include "local_poincare4_projection.hpp"
#include "local_poincare4_projection_spec.hpp"

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
        bool src_specialized,
        bool dst_specialized)
            : m_vector_field(vector_field)
            , m_constraint(constraint)
            , m_order(order)
            , m_src_specialized(src_specialized)
            , m_dst_specialized(dst_specialized)
            , m_src_coordsys(specialize_coordsys(src_coordsys, src_specialized))
            , m_dst_coordsys(specialize_coordsys(dst_coordsys, dst_specialized))
    {
        assert_with_exception(m_vector_field.dimension() == 4);
        assert_with_exception(m_vector_field.imageDimension() == 4);
        assert_with_exception(m_constraint.dimension() == 4);
        assert_with_exception(m_constraint.imageDimension() == 1);
        assert_with_exception(m_src_coordsys.get_origin().dimension() == 4);
        assert_with_exception(m_dst_coordsys.get_origin().dimension() == 4);

        auto is_u_v_pu_zero = [](VectorType origin) -> bool
        {
            return origin[0] == 0.0 && origin[1] == 0.0 && origin[2] == 0.0;
        };

        assert_with_exception(src_specialized == is_u_v_pu_zero(m_src_coordsys.get_origin()));
        assert_with_exception(dst_specialized == is_u_v_pu_zero(m_dst_coordsys.get_origin()));
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
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Specialize local coordinate system for psi_0
    //! @details Adjust Poincare section normal vector so that the section is { v == 0 }.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static CapdUtils::LocalCoordinateSystem<MapT> specialize_coordsys(
        const CapdUtils::LocalCoordinateSystem<MapT>& coordsys,
        bool specialize)
    {
        if (specialize)
        {
            MatrixType new_direction_matrix = coordsys.get_directions_matrix();
            new_direction_matrix(1, 3) = 0.0;
            new_direction_matrix(2, 3) = 1.0;
            new_direction_matrix(3, 3) = 0.0;
            new_direction_matrix(4, 3) = 0.0;
            CapdUtils::LocalCoordinateSystem<MapT> ret
            {
                coordsys.get_origin(),
                new_direction_matrix
            };
            return ret;
        }
        
        return coordsys;
    }


    MapT& m_vector_field;
    MapT& m_constraint;

    unsigned m_order;

    const bool m_src_specialized;
    const bool m_dst_specialized;
    
    const CapdUtils::LocalCoordinateSystem<MapT> m_src_coordsys;
    const CapdUtils::LocalCoordinateSystem<MapT> m_dst_coordsys;

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
            return m_src_specialized ?
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

    using LocalPoincare4_Projection_BaseType = LocalPoincare4_ProjectionBase<MapT>;
    using LocalPoincare4_Projection_BaseTypePtr = std::unique_ptr<LocalPoincare4_Projection_BaseType>;
    using LocalPoincare4_Projection_Type = LocalPoincare4_Projection<MapT>;
    using LocalPoincare4_Projection_SpecType = LocalPoincare4_Projection_Spec<MapT>;

    LocalPoincare4_Projection_BaseTypePtr m_projection_to_2_ptr
    {
        [this]() -> LocalPoincare4_Projection_BaseTypePtr
        {
            return m_dst_specialized ?
                LocalPoincare4_Projection_BaseTypePtr(std::make_unique<LocalPoincare4_Projection_SpecType>( std::ref(m_dst_coordsys) ) ) : 
                LocalPoincare4_Projection_BaseTypePtr(std::make_unique<LocalPoincare4_Projection_Type>( std::ref(m_dst_coordsys) ) );
        }()
    };

    LocalPoincare4_Projection_BaseType& m_projection_to_2
    {
        *m_projection_to_2_ptr
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
