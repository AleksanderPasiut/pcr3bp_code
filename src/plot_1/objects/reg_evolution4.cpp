///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "reg_evolution4.hpp"

#include <pcr3bp_basic/regularized_system.hpp>

namespace Pcr3bpProof
{

using MapT = RMap;

static MapT create_map(Pcr3bp::SetupParameters<MapT> setup, bool positive)
{
    return positive ?
        Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(2, setup, false) :
        Pcr3bp::RegularizedSystem<MapT>::createNegativeVectorField(2, setup, false);
}

static Leo::LinearInterpolationNodesList<double, double>
compute_solution_interpolation_nodes(CapdUtils::SolutionCurve<RMap> const & solution, Leo::Ruler<double> const & ruler)
{
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;

    Leo::LinearInterpolationNodesList<double, double> length_time_nodes = {};

    ScalarType length = 0;
    for (size_t n = 0; n < ruler.count() - 1; ++n)
    {
        const ScalarType time = ruler.minimum() + n * ruler.step();

        const VectorType dU = solution.timeDerivative(time);
        const ScalarType dL = sqrt(1 + capd::vectalg::scalarProduct(dU, dU));

        length_time_nodes.emplace_back( length, time );

        length += dL;
    }

    length_time_nodes.emplace_back( length, ruler.maximum() );

    return length_time_nodes;
}

RegEvolution4::RegEvolution4(
    Lyra::Core3dObjects& core_objects_ref,
    const Manifold4_Transformation & transformation_ref,
    const Params& params)
        : m_transformation_ref(transformation_ref)
        , m_params(params)
        , m_ruler( Leo::Ruler<ScalarType>(0.0, params.t, params.point_count, 1) )
        , m_renderable(
            core_objects_ref,
            Leo::RulerSet<1>({ m_ruler }))
{
    MapT map = create_map(params.setup, params.positive);
    CapdUtils::TimemapWrapper<MapT> timemap( map, 0.0, 20 );

    timemap.set_time(m_params.t);
    timemap(m_params.initial_point, m_solution);

    SolutionInterpolationNodes length_time_nodes = compute_solution_interpolation_nodes(m_solution, m_ruler);

    m_length_time_interpolation = std::make_unique<SolutionInterpolation>( length_time_nodes );

    refresh();
}

void RegEvolution4::refresh()
{
    auto func = [this](const std::array<double, 1>& in) -> std::array<float, 4>
    {
        const ScalarType arg = std::min(in[0] / m_solution.getRightDomain(), 1.0) * m_length_time_interpolation->get_max_argument();
        const ScalarType t = m_length_time_interpolation->evaluate( arg ); 

        const VectorType U = m_solution(t);

        std::array<double, 4> tmp = { U[2], U[3], U[0], U[1] };
        Lyra::Point4d out = m_transformation_ref.func()(tmp);
        return out;
    };

    m_renderable.fill( func, m_params.thickness );
}

}
