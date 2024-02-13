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

template<typename MapT>
class LocalPoincare4_ProjectionBase : public CapdUtils::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    virtual VectorType operator() (const VectorType& vec) override = 0;

    virtual VectorType operator() (const VectorType& vec, MatrixType& der) override = 0;

    virtual unsigned dimension() const override = 0;

    virtual unsigned imageDimension() const override = 0;
};

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

    MatrixType m_dst_coordsys_direction_matrix_inverse
    {
        CapdUtils::gaussInverseMatrix<MapT>( m_dst_coordsys_4_dim.get_directions_matrix() )
    };

    CapdUtils::AffineMap<MapT> m_affine_map
    {
        m_dst_coordsys_4_dim.get_origin(),
        m_dst_coordsys_4_dim.get_directions_matrix()
        // m_dst_coordsys_direction_matrix_inverse
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

    static CapdUtils::LocalCoordinateSystem<MapT> specialize_coordsys(const CapdUtils::LocalCoordinateSystem<MapT>& coordsys, bool specialize)
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

        print_var(src_specialized);
        print_var(dst_specialized);
        print_var(m_src_coordsys.get_origin());
        print_var(m_dst_coordsys.get_origin());
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
