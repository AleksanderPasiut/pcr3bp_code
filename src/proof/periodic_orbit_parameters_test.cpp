///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"
#include "tools/print_bootstrap.hpp"
#include "tools/floating_info.hpp"

#include <capd_utils/map_base.hpp>
#include <capd_utils/pne_map.hpp>
#include <capd_utils/identity_map.hpp>
#include <capd_utils/readable_interval.hpp>
#define CARINA_LOG
#include <capd_utils/newton_method/newton_method.hpp>
#undef CARINA_LOG
#include <capd_utils/poincare_wrapper.hpp>
#include <capd_utils/timemap_wrapper.hpp>
#include <capd_utils/enp_map.hpp>
#include <capd_utils/extension_map.hpp>
#include <capd_utils/composite_map.hpp>
#include <capd_utils/parallel_shooting/parallel_shooting.hpp>
#include <capd_utils/parallel_shooting/parallel_shooting_init.hpp>

#include "pcr3bp_reg_basic_objects.hpp"

#include "tools/variable_printer.hpp"

namespace Pcr3bpProof
{

template<typename MapT>
class LyapunovOrbitRegCollisionSetup
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LyapunovOrbitRegCollisionSetup(ScalarType intermediate_time, bool full_test = true) : m_intermediate_time(intermediate_time)
    {
        std::cout << "Intermediate time: " << intermediate_time << '\n';

        const VectorType h0 = VectorType{ -0.711059 };
        const VectorType v = m_parallel_shooting_init(h0);
        const VectorType v_optimized = VectorType{ v[1], v[2], v[3], v[4], v[0] };

        CapdUtils::NewtonMethod newton(m_parallel_shooting_optimized, v_optimized, 100);
        EXPECT_TRUE(newton.is_successful());

        const VectorType root = newton.get_root();

        const VectorType intermediate_point = { root[0], root[1], root[2], root[3] };
        VectorType image_point = CapdUtils::Extract<MapT>::get_vector(m_poincare(root), 0, 4);

        image_point[1] = 0.0; // v-component is 0.0 from the definition of the poincare section
        image_point[2] = 0.0; // pu-component is 0.0 from the Newton operator definition

        if (full_test)
        {
            EXPECT_EQ(intermediate_point, m_basic_objects.m_parameters.get_intermediate_point());
            EXPECT_EQ(image_point, m_basic_objects.m_parameters.get_image_point());

            {
                const ScalarType lyapunov_orbit_period_diff = capd::abs( (intermediate_time + m_poincare.get_last_evaluation_return_time())*2 - m_basic_objects.m_lyapunov_orbit_period );
                EXPECT_LT(lyapunov_orbit_period_diff, 1.3e-13);
            }

            const std::string suffix = std::is_same<IMap, MapT>::value ? ".txt" : ".approx.txt";

            CapdUtils::VariablePrinter<MapT>::print(
                "periodic_orbit_initial_point" + suffix,
                "Initial point on periodic orbit, (u,v,pu,pv) components",
                m_initial_point);

            CapdUtils::VariablePrinter<MapT>::print(
                "periodic_orbit_intermediate_point" + suffix,
                "Intermediate point on periodic orbit, (u,v,pu,pv) components",
                intermediate_point);

            CapdUtils::VariablePrinter<MapT>::print(
                "periodic_orbit_image_point" + suffix,
                "Image point on periodic orbit, (u,v,pu,pv) components",
                image_point);

            const ScalarType energy = root[4];
            CapdUtils::VariablePrinter<MapT>::print(
                "collision_energy" + suffix,
                "Energy of Lyapunov orbit that passes through regularized collision",
                energy);
        }
    }

private:
    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    ScalarType const m_intermediate_time;
    MapT m_vf_reg { Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField(2, m_basic_objects.m_setup, false) };

    CapdUtils::CoordinateSection<MapT> m_section { m_vf_reg.dimension(), 1, 0 };
    CapdUtils::PoincareWrapper<MapT, decltype(m_section)> m_poincare { m_vf_reg, m_basic_objects.m_order, m_section };
    CapdUtils::TimemapWrapper<MapT> m_timemap { m_vf_reg, m_intermediate_time, m_basic_objects.m_order };

    const VectorType m_initial_point { m_basic_objects.m_parameters.get_initial_point() };
    const VectorType m_h_to_5_extend_vector { CapdUtils::Concat<MapT>::concat_vectors({ m_initial_point, VectorType{ 0.0 } }) };
    MapT m_h_to_5_extend { CapdUtils::ExtensionMap<MapT>::create(m_h_to_5_extend_vector, { -1, -1, -1, -1, 0 }) };
    CapdUtils::CompositeMap<MapT, MapT, decltype(m_timemap)&> m_timemap_1_to_5
    {
        m_h_to_5_extend,
        std::ref(m_timemap)
    };

    MapT m_5_to_pu_projection { CapdUtils::ProjectionMap<MapT>::create(5, { 2 }) };
    CapdUtils::CompositeMap<MapT, decltype(m_poincare)&, MapT> m_poincare_5_to_1
    {
        std::ref(m_poincare),
        m_5_to_pu_projection
    };

    CapdUtils::ParallelShootingInit<MapT, decltype(m_timemap_1_to_5)&, decltype(m_poincare_5_to_1)&> m_parallel_shooting_init
    {
        std::ref(m_timemap_1_to_5),
        std::ref(m_poincare_5_to_1)
    };

    CapdUtils::ParallelShooting<MapT, decltype(m_timemap_1_to_5)&, decltype(m_poincare_5_to_1)&> m_parallel_shooting
    {
        std::ref(m_timemap_1_to_5),
        std::ref(m_poincare_5_to_1)
    };

    CapdUtils::ENP<MapT, decltype(m_parallel_shooting)&> m_parallel_shooting_optimized
    {
        VectorType{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 4, 0, 1, 2, 3, 4 },
        { 0, 1, 2, 3, 5 },
        std::ref(m_parallel_shooting)
    };
};

}

TEST(Pcr3bp_intermediate, periodic_orbit_parameters_test_nonrigorous)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();
    LyapunovOrbitRegCollisionSetup<RMap> setup( 58696.0 / 65536 );
}

TEST(Pcr3bp_intermediate, periodic_orbit_parameters_test_rigorous)
{
    using namespace Pcr3bpProof;

    LyapunovOrbitRegCollisionSetup<IMap> setup( 58696.0 / 65536 );
}

TEST(Extended_Pcr3bp_intermediate, periodic_orbit_parameters_test_rigorous_2)
{
    using namespace Pcr3bpProof;

    LyapunovOrbitRegCollisionSetup<IMap> setup( 7.0 / 8, false );
}

TEST(Exteded_Pcr3bp_intermediate, periodic_orbit_parameters_test_rigorous_3)
{
    using namespace Pcr3bpProof;

    LyapunovOrbitRegCollisionSetup<IMap> setup( 1.0 / 2, false );
}
