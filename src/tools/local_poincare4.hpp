///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"

#include <carina/poincare_wrapper.hpp>
#include <carina/timemap_wrapper.hpp>
#include <carina/extension_map.hpp>
#include <carina/projection_map.hpp>
#include <carina/constrained_function.hpp>

#include "id_with_constraint.hpp"
#include "affine_poincare_map.hpp"

#include "local_poincare4_constraint.hpp"

namespace Ursa
{

template<typename MapT>
class LocalPoincare4 : public Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LocalPoincare4(
        MapT& vector_field,
        MapT& constraint,
        unsigned order,
        const Carina::LocalCoordinateSystem<MapT>& src_coordsys,
        const Carina::LocalCoordinateSystem<MapT>& dst_coordsys)
            : m_vector_field(vector_field)
            , m_constraint(constraint)
            , m_order(order)
            , m_src_coordsys(src_coordsys)
            , m_dst_coordsys(dst_coordsys)
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

    void operator() (const VectorType& vec, ScalarType time, Carina::SolutionCurve<MapT>& solution_curve)
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
    
    const Carina::LocalCoordinateSystem<MapT> m_src_coordsys;
    const Carina::LocalCoordinateSystem<MapT> m_dst_coordsys;

    Carina::AffinePoincareMap<MapT> m_affine_poincare
    {
        m_vector_field,
        m_order,
        m_src_coordsys,
        m_dst_coordsys,
    };

    LocalPoincare4_Constraint<MapT> m_extension_to_4
    {
        std::ref(m_constraint),
        std::ref(m_src_coordsys)
    };

    MapT m_projection_to_2
    {
        Carina::ProjectionMap<MapT>::create( 4, { 0, 1 } )
    };

    Carina::CompositeMap<MapT,
        decltype(m_extension_to_4)&,
        decltype(m_affine_poincare)&,
        decltype(m_projection_to_2)&> m_affine_poincare_2
    {
        std::ref(m_extension_to_4),
        std::ref(m_affine_poincare),
        std::ref(m_projection_to_2)
    };

    Carina::AffineMap<MapT> m_affine_src
    {
        m_src_coordsys
    };

    Carina::TimemapWrapper<MapT> m_timemap
    {
        m_vector_field,
        0.0,
        m_order
    };
};

}
