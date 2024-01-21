///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/local_poincare4.hpp"
#include "tools/gain_map.hpp"

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Implementation of the local poincare map with scaling factor included
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class ScaledLocalPoincare4_Map : public CapdUtils::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    ScaledLocalPoincare4_Map(
        MapT& vector_field,
        MapT& constraint,
        unsigned order,
        const CapdUtils::LocalCoordinateSystem<MapT>& src_coordsys,
        const CapdUtils::LocalCoordinateSystem<MapT>& dst_coordsys,
        ScalarType input_gain,
        bool src_specialized = false,
        bool dst_specialized = false)
            : m_local_poincare4(
                vector_field,
                constraint,
                order,
                src_coordsys,
                dst_coordsys,
                src_specialized,
                dst_specialized)
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
    void operator() (const VectorType& vec, ScalarType time, CapdUtils::SolutionCurve<MapT>& solution_curve)
    {
        const ScalarType k = m_input_gain.get_gain();
        return m_local_poincare4(vec * k, time, solution_curve);
    }

private:
    LocalPoincare4<MapT> m_local_poincare4;
    CapdUtils::GainMap<MapT> m_input_gain;
    CapdUtils::GainMap<MapT> m_output_gain;

    CapdUtils::CompositeMap<
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
