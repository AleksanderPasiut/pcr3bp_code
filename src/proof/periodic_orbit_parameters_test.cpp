///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"
#include "tools/print_bootstrap.hpp"
#include "tools/floating_info.hpp"

#include <carina/map_base.hpp>
#include <carina/pne_map.hpp>
#include <carina/identity_map.hpp>
#include <carina/readable_interval.hpp>
#define CARINA_LOG
#include <carina/newton_method/newton_method.hpp>
#undef CARINA_LOG
#include <carina/poincare_wrapper.hpp>
#include <carina/timemap_wrapper.hpp>
#include <carina/enp_map.hpp>
#include <carina/extension_map.hpp>
#include <carina/composite_map.hpp>
#include <carina/parallel_shooting/parallel_shooting.hpp>
#include <carina/parallel_shooting/parallel_shooting_init.hpp>

#include "pcr3bp_reg_basic_objects.hpp"

#include "tools/variable_printer.hpp"

namespace Ursa
{

template<typename MapT>
class LyapunovOrbitRegCollisionSetup
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LyapunovOrbitRegCollisionSetup()
    {
        const VectorType h0 = VectorType{ -0.711059 };
        const VectorType v = m_parallel_shooting_init(h0);
        const VectorType v_optimized = VectorType{ v[1], v[2], v[3], v[4], v[0] };

        Carina::NewtonMethod newton(m_parallel_shooting_optimized, v_optimized, 100);
        EXPECT_TRUE(newton.is_successful());

        const VectorType root = newton.get_root();

        const VectorType intermediate_point = { root[0], root[1], root[2], root[3] };
        VectorType image_point = Carina::Extract<MapT>::get_vector(m_poincare(root), 0, 4);

        image_point[1] = 0.0; // v-component is 0.0 from the definition of the poincare section
        image_point[2] = 0.0; // pu-component is 0.0 from the Newton operator definition

        EXPECT_EQ(intermediate_point, m_basic_objects.m_parameters.get_intermediate_point());
        EXPECT_EQ(image_point, m_basic_objects.m_parameters.get_image_point());

        const std::string suffix = std::is_same<IMap, MapT>::value ? ".txt" : ".approx.txt";

        Carina::VariablePrinter<MapT>::print(
            "periodic_orbit_initial_point" + suffix,
            "Initial point on periodic orbit, (u,v,pu,pv) components",
            m_initial_point);

        Carina::VariablePrinter<MapT>::print(
            "periodic_orbit_intermediate_point" + suffix,
            "Intermediate point on periodic orbit, (u,v,pu,pv) components",
            intermediate_point);

        Carina::VariablePrinter<MapT>::print(
            "periodic_orbit_image_point" + suffix,
            "Image point on periodic orbit, (u,v,pu,pv) components",
            image_point);

        const ScalarType energy = root[4];
        Carina::VariablePrinter<MapT>::print(
            "collision_energy" + suffix,
            "Energy of Lyapunov orbit that passes through regularized collision",
            energy);
    }

private:
    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    ScalarType const m_t0 { 58696.0 / 65536 };
    MapT m_vf_reg { Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(2, m_basic_objects.m_setup, false) };

    Carina::CoordinateSection<MapT> m_section { m_vf_reg.dimension(), 1, 0 };
    Carina::PoincareWrapper<MapT, decltype(m_section)> m_poincare { m_vf_reg, m_basic_objects.m_order, m_section };
    Carina::TimemapWrapper<MapT> m_timemap { m_vf_reg, m_t0, m_basic_objects.m_order };

    const VectorType m_initial_point { m_basic_objects.m_parameters.get_initial_point() };
    const VectorType m_h_to_5_extend_vector { Carina::Concat<MapT>::concat_vectors({ m_initial_point, VectorType{ 0.0 } }) };
    MapT m_h_to_5_extend { Carina::ExtensionMap<MapT>::create(m_h_to_5_extend_vector, { -1, -1, -1, -1, 0 }) };
    Carina::CompositeMap<MapT, MapT, decltype(m_timemap)&> m_timemap_1_to_5
    {
        m_h_to_5_extend,
        std::ref(m_timemap)
    };

    MapT m_5_to_pu_projection { Carina::ProjectionMap<MapT>::create(5, { 2 }) };
    Carina::CompositeMap<MapT, decltype(m_poincare)&, MapT> m_poincare_5_to_1
    {
        std::ref(m_poincare),
        m_5_to_pu_projection
    };

    Carina::ParallelShootingInit<MapT, decltype(m_timemap_1_to_5)&, decltype(m_poincare_5_to_1)&> m_parallel_shooting_init
    {
        std::ref(m_timemap_1_to_5),
        std::ref(m_poincare_5_to_1)
    };

    Carina::ParallelShooting<MapT, decltype(m_timemap_1_to_5)&, decltype(m_poincare_5_to_1)&> m_parallel_shooting
    {
        std::ref(m_timemap_1_to_5),
        std::ref(m_poincare_5_to_1)
    };

    Carina::ENP<MapT, decltype(m_parallel_shooting)&> m_parallel_shooting_optimized
    {
        VectorType{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 4, 0, 1, 2, 3, 4 },
        { 0, 1, 2, 3, 5 },
        std::ref(m_parallel_shooting)
    };
};

}

TEST(Pcr3bp_intermediate, periodic_orbit_parameters_test_rigorous)
{
    using namespace Ursa;

    LyapunovOrbitRegCollisionSetup<IMap> setup {};
}
