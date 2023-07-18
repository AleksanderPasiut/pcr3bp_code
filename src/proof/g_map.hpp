///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/local_poincare4.hpp"
#include "tools/gig_map.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Implementation of the `g` map
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class G_Map : public Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    G_Map(
        MapT& vector_field,
        MapT& constraint,
        unsigned order,
        const Carina::LocalCoordinateSystem<MapT>& src_coordsys,
        const Carina::LocalCoordinateSystem<MapT>& dst_coordsys,
        ScalarType input_gain)
            : m_local_poincare4(
                vector_field,
                constraint,
                order,
                src_coordsys,
                dst_coordsys)
            , m_gig(m_local_poincare4, input_gain)
    {}

    VectorType operator() (const VectorType& vec) override
    {
        return m_gig(vec);
    }

    VectorType operator() (const VectorType& vec, MatrixType& der) override
    {
        return m_gig(vec, der);
    }

    unsigned dimension() const override
    {
        return m_gig.dimension();
    }

    unsigned imageDimension() const override
    {
        return m_gig.imageDimension();
    }

    ScalarType get_last_evaluation_return_time() const
    {
        return m_local_poincare4.get_last_evaluation_return_time();
    }

    void operator() (const VectorType& vec, ScalarType time, Carina::SolutionCurve<MapT>& solution_curve)
    {
        return m_local_poincare4(vec * get_input_gain(), time, solution_curve);
    }

    ScalarType get_input_gain()
    {
        return m_gig.get_input_gain();
    }

private:
    LocalPoincare4<MapT> m_local_poincare4;
    Carina::GigMap<MapT, decltype(m_local_poincare4)> m_gig;
};

}
