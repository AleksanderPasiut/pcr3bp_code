///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>

#include "tools/plotting/solution_curve_interpolation.hpp"
#include "capd_renderable.hpp"

#include "ppv_manifolds_maps.hpp"

namespace Ursa
{

template<typename MapT>
class PpvManifold
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using InterpolationNodes = Leo::CubicInterpolationNodesList<ScalarType, VectorType>;

    PpvManifold(const InterpolationNodes& nodes)
        : m_nodes(nodes)
        , m_interpolation(m_nodes)
    {}

    VectorType operator() (const VectorType& t) const
    {
        const ScalarType arg = std::min(t[0], m_interpolation.get_max_argument());
        const VectorType ret = m_interpolation.evaluate(arg);
        return ret;
    }

    ScalarType get_min_argument() const
    {
        return m_interpolation.get_min_argument();
    }

    ScalarType get_max_argument() const
    {
        return m_interpolation.get_max_argument();
    }

    const InterpolationNodes& get_nodes() const noexcept
    {
        return m_nodes;
    }

    std::list<VectorType> find_intersections(const PpvManifold<MapT>& arg) const
    {
        std::list<VectorType> ret;

        auto l = m_interpolation.find_intersections(arg.m_interpolation, 5.0e-5, 0.0);

        for (auto & x : l)
        {
            ret.emplace_back(x.point);
        }

        return ret;
    }

private:
    InterpolationNodes m_nodes;
    Leo::CubicInterpolationStandard<ScalarType, VectorType> m_interpolation;
};

struct PpvManifoldInitParam
{
    RVector origin;
    RVector direction;
    size_t points;
    double min_multiplier;
    double max_multiplier;
};

template<typename PoincareT>
class PpvManifoldFactory
{
public:
    PpvManifoldFactory(const Pcr3bpSetupValues<RMap>& setup)
        : m_rhez_u_alpha(setup)
        , m_poincare(setup)
    {}

    PpvManifold<RMap> init(const PpvManifoldInitParam& param)
    {
        PpvManifold<RMap>::InterpolationNodes nodes {};

        const double a = (param.min_multiplier - param.max_multiplier)/2.0;
        const double b = (param.max_multiplier + param.min_multiplier)/2.0;

        for (size_t i = 0; i < param.points; ++i)
        {
            const double dx = a*cos(M_PI*double(2*i+1)/(2*param.points)) + b;

            const RVector ptr = param.origin + param.direction * dx;

            RMatrix mat(5, 2);
            const RVector qtr = m_rhez_u_alpha(ptr, mat);
            const RVector der = mat * param.direction;

            nodes.emplace_back( dx, qtr, der );
        }

        return PpvManifold<RMap>(nodes);
    }

    PpvManifold<RMap> evolve(const PpvManifold<RMap>& prev)
    {
        PpvManifold<RMap>::InterpolationNodes nodes {};

        for (const auto& node : prev.get_nodes())
        {
            RMatrix mat(5, 5);
            const RVector qtr = m_poincare(node.value, mat);
            const RVector der = mat * node.derivative;

            nodes.emplace_back( node.argument, qtr, der );
        }

        return PpvManifold<RMap>(nodes);
    }

    // PpvManifold<RMap> reduce(const PpvManifold<RMap>& prev, double min_node_dist)
    // {
    //     PpvManifold<RMap>::InterpolationNodes nodes {};
    //
    //     Carina::MaxNorm<RMap> norm;
    //
    //     const Leo::CubicInterpolationNode<double, RVector>* last_node = nullptr;
    //     for (const auto& node : prev.get_nodes())
    //     {
    //         const double d = last_node ? norm(node.value - last_node->value) : 100;
    //         if ((last_node == nullptr) || (d > min_node_dist))
    //         {
    //             nodes.emplace_back( node.argument, node.value, node.derivative );
    //             last_node = &node;
    //         }
    //     }
    //
    //     if (nodes.size() <= 2)
    //     {
    //         return prev;
    //     }
    //     else
    //     {
    //         return PpvManifold<RMap>(nodes);
    //     }
    // }

private:
    RhezUAlpha_Render2<RMap> m_rhez_u_alpha;
    PoincareT m_poincare;
};

}
