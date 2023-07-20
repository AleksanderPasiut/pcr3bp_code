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
            , m_input_gain(input_gain, m_local_poincare4.dimension())
            , m_output_gain(ScalarType(1.0) / input_gain, m_local_poincare4.imageDimension())
    {}

    VectorType operator() (const VectorType& vec) override
    {
        return m_composite(vec);
    }

    VectorType operator() (const VectorType& vec, MatrixType& der) override
    {
        return m_composite(vec, der);
    }

    unsigned dimension() const override
    {
        return m_composite.dimension();
    }

    unsigned imageDimension() const override
    {
        return m_composite.imageDimension();
    }

    ScalarType get_last_evaluation_return_time() const
    {
        return m_local_poincare4.get_last_evaluation_return_time();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Get Poincare time and solution curve of the underlying Poincare map
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void operator() (const VectorType& vec, ScalarType time, Carina::SolutionCurve<MapT>& solution_curve)
    {
        const ScalarType k = m_input_gain.get_gain();
        return m_local_poincare4(vec * k, time, solution_curve);
    }

    // ScalarType get_input_gain()
    // {
    //     return m_input_gain.get_gain();
    // }

private:
    LocalPoincare4<MapT> m_local_poincare4;
    Carina::GainMap<MapT> m_input_gain;
    Carina::GainMap<MapT> m_output_gain;

    Carina::CompositeMap<
        MapT,
        decltype(m_input_gain)&,
        decltype(m_local_poincare4)&,
        decltype(m_output_gain)&> m_composite
    {
        std::ref(m_input_gain),
        std::ref(m_local_poincare4),
        std::ref(m_output_gain)
    };
};

}
