///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/map_base.hpp>
#include <carina/extract.hpp>
#include <carina/affine_map.hpp>
#include <carina/local_poincare_wrapper.hpp>

namespace Carina
{

template<typename MapT>
class AffinePoincareMap : public Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    AffinePoincareMap(
        MapT& vector_field,
        unsigned order,
        LocalCoordinateSystem<MapT> src_coordsys,
        LocalCoordinateSystem<MapT> dst_coordsys)
            : m_src_linear_map( src_coordsys )
            , m_dst_section( gen_section(dst_coordsys) )
            , m_poincare( vector_field, order, m_dst_section, src_coordsys, dst_coordsys )
    {
        assert_with_exception(vector_field.dimension() == src_coordsys.get_origin().dimension());
        assert_with_exception(vector_field.dimension() == dst_coordsys.get_origin().dimension());
        assert_with_exception(vector_field.dimension() >= 3);
    }

    VectorType operator() (const VectorType& vec) override
    {
        return m_poincare(vec);
    }

    VectorType operator() (const VectorType& vec, MatrixType& der) override
    {
        return m_poincare(vec, der);
    }

    unsigned dimension() const override
    {
        return m_poincare.dimension();
    }

    unsigned imageDimension() const override
    {
        return m_poincare.imageDimension();
    }

    ScalarType get_last_evaluation_return_time() const
    {
        return m_poincare.get_last_evaluation_return_time();
    }

private:
    static AffineSection<MapT> gen_section( const LocalCoordinateSystem<MapT>& coordsys )
    {
        const VectorType vector_field_dir = Extract<MapT>::get_vvector(coordsys.get_directions_matrix(), 3);
        return AffineSection<MapT>(coordsys.get_origin(), vector_field_dir);            
    }

    AffineMap<MapT> m_src_linear_map;
    AffineSection<MapT> m_dst_section;

    LocalPoincareWrapper<MapT, AffineSection<MapT>> m_poincare;
};

}
