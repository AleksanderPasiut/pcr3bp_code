///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"
#include <capd_utils/capd/solution_curve.hpp>
#include <leo/interpolation/interpolation.hpp>

namespace Ursa
{

template<typename MapT>
class SolutionCurveInterpolation
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    SolutionCurveInterpolation(const Carina::SolutionCurve<MapT>& solution, size_t nodes_count)
        : m_interpolation(get_nodes(solution, nodes_count))
    {}

    VectorType operator() (const VectorType& vec)
    {
        return m_interpolation.evaluate(vec[0]);
    }

private:
    Leo::CubicInterpolationNodesList<ScalarType, VectorType> get_nodes(
        const Carina::SolutionCurve<MapT>& solution,
        size_t nodes_count)
    {
        Leo::CubicInterpolationNodesList<ScalarType, VectorType> nodes;

        const Real dt = (solution.getRightDomain() - solution.getLeftDomain()) / (nodes_count - 1);

        Real t = 0;
        for (size_t i = 1; i < nodes_count; ++i, t += dt)
        {
            const VectorType val = solution(t);
            const VectorType der = solution.timeDerivative(t);

            nodes.emplace_back(ScalarType(t), val, der);
        }

        {
            const Real t = solution.getRightDomain();
            const VectorType val = solution(t);
            const VectorType der = solution.timeDerivative(t);
            nodes.emplace_back(ScalarType(t), val, der);
        }

        return nodes;
    }

    Leo::CubicInterpolationStandard<ScalarType, VectorType> m_interpolation;
};

}
