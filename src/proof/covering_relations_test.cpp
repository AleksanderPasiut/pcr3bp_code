///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"

#include "tools/coordsys_utilities.hpp"
#include "tools/solution_curve_with_condition_check.hpp"
#include "tools/auxiliary_functions.hpp"

#include "covering_relations_setup.hpp"
#include "covering_relation_checker.hpp"

#include "g_map.hpp"

namespace Pcr3bpProof
{

template<typename MapT>
class CoveringRelationsTest
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = CapdUtils::LocalCoordinateSystem<MapT>;

    CoveringRelationsTest(const CoveringRelationsSetup& setup)
        : m_periodic_orbit_coordsys(setup.get_periodic_orbit_coordsys())
        , m_homoclinic_orbit_coordsys(setup.get_homoclinic_orbit_coordsys())
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check covering relations along homoclinic orbit
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void check_homoclinic_coverings()
    {
        for (size_t i = 1; i < m_homoclinic_orbit_coordsys.size(); ++i)
        {
            const size_t src_idx = i-1;
            const size_t dst_idx = i;
            std::cout << "homoclinic orbit covering " << src_idx << " => " << dst_idx << '\n';

            const CapdUtils::LocalCoordinateSystem<MapT> coordsys_src = m_homoclinic_orbit_coordsys.at(src_idx);
            const CapdUtils::LocalCoordinateSystem<MapT> coordsys_dst = m_homoclinic_orbit_coordsys.at(dst_idx);

            const ScalarType time_span = check_covering_relation_forward(coordsys_src, coordsys_dst);
            simple_collision_avoidance_check(coordsys_src, coordsys_dst, time_span);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check covering relations along periodic orbit
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void check_periodic_coverings()
    {
        {
            std::cout << "periodic orbit covering 0 => 1\n";

            const CapdUtils::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(0);
            const CapdUtils::LocalCoordinateSystem<MapT> coordsys_dst = m_periodic_orbit_coordsys.at(1);
            check_covering_relation_forward(coordsys_src, coordsys_dst);
        }

        {
            std::cout << "periodic orbit covering 1 => 2\n";

            const CapdUtils::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(1);
            const CapdUtils::LocalCoordinateSystem<MapT> coordsys_dst = m_periodic_orbit_coordsys.at(2);
            const ScalarType time_span = check_covering_relation_forward(coordsys_src, coordsys_dst);
            simple_collision_avoidance_check(coordsys_src, coordsys_dst, time_span);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check covering relations between homoclinic and periodic orbits
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void check_jump_coverings()
    {
        {
            const CapdUtils::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(3);
            const CapdUtils::LocalCoordinateSystem<MapT> coordsys_dst = *( m_homoclinic_orbit_coordsys.begin() );

            std::cout << "periodic (3) <= first homoclinic covering\n";

            const ScalarType time_span = check_covering_relation_backward(coordsys_src, coordsys_dst);
            simple_collision_avoidance_check(coordsys_src, coordsys_dst, time_span);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check parallelogram covering relations
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void check_parallelogram_coverings()
    {
        parallelogram_covering_check();
        parallelogram_covering_endings_check();
        collision_manifold_derivative();
    }

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check forward covering relation between given local coordinate systems
    //! @return Time interval of underlying evolved trajectory
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ScalarType check_covering_relation_forward(
        CapdUtils::LocalCoordinateSystem<MapT> coordsys_src,
        CapdUtils::LocalCoordinateSystem<MapT> coordsys_dst)
    {
        G_Map<MapT> f
        {
            std::ref(m_basic_objects.m_vf_reg_pos2),
            std::ref(m_basic_objects.m_hamiltonian_reg2),
            m_basic_objects.m_order,
            coordsys_src,
            coordsys_dst,
            m_gain_factor
        };

        CoveringRelationCheck cr { f };

        const ScalarType time_span = f.get_last_evaluation_return_time();

        EXPECT_TRUE(cr.contraction_condition());
        EXPECT_TRUE(cr.expansion_condition());

        // check that image is properly covered by its coordinate system
        LocalPoincare4_Constraint<MapT> extension_to_4_dst
        {
            std::ref(m_basic_objects.m_hamiltonian_reg2),
            std::ref(coordsys_dst)
        };

        extension_to_4_dst(cr.get_img() * m_gain_factor);
        // print_var( extension_to_4_dst(img) );

        return time_span;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check backward covering relation between given local coordinate systems
    //! @return Return interval of time of evolved solution
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ScalarType check_covering_relation_backward(
        CapdUtils::LocalCoordinateSystem<MapT> coordsys_src,
        CapdUtils::LocalCoordinateSystem<MapT> coordsys_dst)
    {
        G_Map<MapT> f
        {
            std::ref(m_basic_objects.m_vf_reg_neg2),
            std::ref(m_basic_objects.m_hamiltonian_reg2),
            m_basic_objects.m_order,
            coordsys_dst,
            coordsys_src,
            m_gain_factor
        };

        MapT J = AuxiliaryFunctions<MapT>::J();
        CapdUtils::CompositeMap<MapT, decltype(J)&, decltype(f)&, decltype(J)&> jfj
        {
            std::ref(J),
            std::ref(f),
            std::ref(J)
        };

        CoveringRelationCheck cr { jfj };

        const ScalarType time_span = f.get_last_evaluation_return_time();

        EXPECT_TRUE(cr.contraction_condition());
        EXPECT_TRUE(cr.expansion_condition());

        // check that image is properly covered by its coordinate system
        LocalPoincare4_Constraint<MapT> extension_to_4_dst
        {
            std::ref(m_basic_objects.m_hamiltonian_reg2),
            std::ref(coordsys_dst)
        };

        extension_to_4_dst(cr.get_img() * m_gain_factor);
        // print_var( extension_to_4_dst(img) );

        return time_span;
    }

    void simple_collision_avoidance_check(
        CapdUtils::LocalCoordinateSystem<MapT> coordsys_src,
        CapdUtils::LocalCoordinateSystem<MapT> coordsys_dst,
        ScalarType time_span)
    {
        LocalPoincare4_Constraint<MapT> extension_to_4_src
        {
            std::ref(m_basic_objects.m_hamiltonian_reg2),
            std::ref(coordsys_src)
        };

        // check that image is properly covered by its coordinate system
        LocalPoincare4_Constraint<MapT> extension_to_4_dst
        {
            std::ref(m_basic_objects.m_hamiltonian_reg2),
            std::ref(coordsys_dst)
        };

        {
            using CapdUtils::Node;
            auto collision_condition_f = [](Node, Node in[], int, Node out[], int, Node param[], int)
            {
                out[0] = in[0];
                out[1] = in[1];
                out[2] = sqr(in[2]) + sqr(in[3]) - 8 * param[0];
            };

            MapT collision_condition(collision_condition_f, 4, 3, 1);
            collision_condition.setParameter(0, m_basic_objects.m_setup.get_mu(2));

            CapdUtils::MaxNorm<MapT> norm {};

            const VectorType expected_collision = m_basic_objects.m_parameters.get_initial_point();
            if ( norm(coordsys_src.get_origin() - expected_collision) < norm(coordsys_dst.get_origin() - expected_collision) )
            {
                G_Map<MapT> f_pos
                {
                    std::ref(m_basic_objects.m_vf_reg_pos2),
                    std::ref(m_basic_objects.m_hamiltonian_reg2),
                    m_basic_objects.m_order,
                    coordsys_src,
                    coordsys_dst,
                    m_gain_factor
                };

                SolutionCurveWithConditionCheck<MapT> solution_curve {};
                f_pos(N, time_span, solution_curve);

                bool const solution_curve_condition = solution_curve.is_condition_never_satisfied( collision_condition );
                EXPECT_TRUE(solution_curve_condition);
            }
            else
            {
                G_Map<MapT> f_neg
                {
                    std::ref(m_basic_objects.m_vf_reg_neg2),
                    std::ref(m_basic_objects.m_hamiltonian_reg2),
                    m_basic_objects.m_order,
                    coordsys_dst,
                    coordsys_src,
                    m_gain_factor
                };

                SolutionCurveWithConditionCheck<MapT> solution_curve {};
                f_neg(N, time_span, solution_curve);

                bool const solution_curve_condition = solution_curve.is_condition_never_satisfied( collision_condition );
                EXPECT_TRUE(solution_curve_condition);
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check parallelogram coverings around fixed point
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void parallelogram_covering_check()
    {
        const VectorType arg = N;
        const ScalarType L = m_basic_objects.m_parallelogram_coverings_parameters.L;

        std::list<MatrixType> der_list {};

        MapT eta = AuxiliaryFunctions<MapT>::eta( L );
        MapT eta_inverse = AuxiliaryFunctions<MapT>::eta( -L );

        for (int i = 0; i < 4; ++i)
        {
            const int first = i;
            const int second = (i+1) % 4;
            
            G_Map<MapT> poincare
            {
                m_basic_objects.m_vf_reg_pos2,
                m_basic_objects.m_hamiltonian_reg2,
                m_basic_objects.m_order,
                m_periodic_orbit_coordsys.at(first),
                m_periodic_orbit_coordsys.at(second),
                m_gain_factor
            };

            CapdUtils::CompositeMap<MapT, MapT&, decltype(poincare)&, MapT&> aligned_poincare
            {
                std::ref(eta),
                std::ref(poincare),
                std::ref(eta_inverse)
            };

            MatrixType der(2,2);
            aligned_poincare(arg, der);
            print_var(der);

            der_list.emplace_back(der);
        }

        auto it = der_list.begin();

        MatrixType der_union = *it;
        for (++it; it != der_list.end(); ++it)
        {
            der_union = capd::vectalg::intervalHull(*it, der_union);
        }

        print_var(der_union);

        const ScalarType alpha = 5.09;
        const ScalarType beta = 0.195;
        const ScalarType rho = 0.197;
        const ScalarType c = 0.0011;

        EXPECT_TRUE( 0 < beta );
        EXPECT_TRUE( beta < rho );
        EXPECT_TRUE( alpha > 2*c + rho );
        EXPECT_TRUE( c + rho < 1 );

        EXPECT_TRUE( der_union(1,1) > alpha );
        EXPECT_TRUE( der_union(1,2) < 0 );
        EXPECT_TRUE( der_union(1,2) > -c );
        EXPECT_TRUE( der_union(2,1) > 0 );
        EXPECT_TRUE( der_union(2,1) < c );
        EXPECT_TRUE( der_union(2,2) < rho );
        EXPECT_TRUE( der_union(2,2) > beta );
    }

    void parallelogram_covering_endings_check()
    {
        const ScalarType L = m_basic_objects.m_parallelogram_coverings_parameters.L;
        const ScalarType b0 = m_basic_objects.m_parallelogram_coverings_parameters.b0;
        const ScalarType a0 = m_basic_objects.m_parallelogram_coverings_parameters.a0;

        EXPECT_TRUE(0 < a0);
        EXPECT_TRUE(a0 < b0);
        EXPECT_TRUE(b0 < 1);

        MapT R_inverse = AuxiliaryFunctions<MapT>::R_Inverse(a0, b0);
        MapT eta_inverse = AuxiliaryFunctions<MapT>::eta( -L );
        MapT J = AuxiliaryFunctions<MapT>::J();

        const CapdUtils::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(3);
        const CapdUtils::LocalCoordinateSystem<MapT> coordsys_dst = *( m_homoclinic_orbit_coordsys.begin() );

        G_Map<MapT> poincare
        {
            m_basic_objects.m_vf_reg_neg2,
            m_basic_objects.m_hamiltonian_reg2,
            m_basic_objects.m_order,
            coordsys_dst,
            coordsys_src,
            m_gain_factor
        };

        CapdUtils::CompositeMap<MapT,
            decltype(J)&,
            decltype(poincare)&,
            decltype(J)&,
            decltype(eta_inverse)&,
            decltype(R_inverse)&> composite
        {
            std::ref(J),
            std::ref(poincare),
            std::ref(J),
            std::ref(eta_inverse),
            std::ref(R_inverse)
        };
        
        CoveringRelationCheck cr { composite };

        EXPECT_TRUE(cr.contraction_condition());
        EXPECT_TRUE(cr.expansion_condition());
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check collision manifold derivative
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void collision_manifold_derivative()
    {
        const ScalarType L = m_basic_objects.m_parallelogram_coverings_parameters.L;

        MapT eta = AuxiliaryFunctions<MapT>::eta( L );

        LocalPoincare4_Constraint<MapT> map_E0
        {
            std::ref(m_basic_objects.m_hamiltonian_reg2),
            std::ref( m_periodic_orbit_coordsys.at(0) )
        };

        CapdUtils::AffineMap<MapT> map_L0 { m_periodic_orbit_coordsys.at(0) };

        using CapdUtils::Node;
        auto collision_condition_f = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            out[0] = sqr(in[2]) + sqr(in[3]) - 8 * param[0];
        };

        MapT map_C(collision_condition_f, 4, 1, 1);
        map_C.setParameter(0, m_basic_objects.m_setup.get_mu(2));

        CapdUtils::CompositeMap<MapT,
            decltype(eta)&,
            decltype(map_E0)&,
            decltype(map_L0)&,
            decltype(map_C)&> map_C2(
                std::ref(eta),
                std::ref(map_E0),
                std::ref(map_L0),
                std::ref(map_C)
            );

        std::cout.precision(20);

        MatrixType der(1, 2);
        map_C2(N * m_gain_factor, der);
        const ScalarType fx = der(1,1);
        const ScalarType fy = der(1,2);

        EXPECT_TRUE( fx < 0 );
        EXPECT_TRUE( fy < 0 );

        std::ofstream fs("collision_manifold_parameters.txt");
        if (fs)
        {
            CapdUtils::VariablePrinter<MapT>::print(fs, "Collision manifold derivative fx", fx);
            CapdUtils::VariablePrinter<MapT>::print(fs, "Collision manifold derivative fy", fy);
            CapdUtils::VariablePrinter<MapT>::print(fs, "Collision manifold derivative -ddy/ddx", -fy / fx);
            CapdUtils::VariablePrinter<MapT>::print(fs, "Collision manifold derivative -ddx/ddy", -fx / fy);
            fs.close();
        }
        else
        {
            throw std::logic_error("Failed to export collision manifold parameters!");
        }
    }


    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    const std::vector<Coordsys> m_periodic_orbit_coordsys;
    const std::vector<Coordsys> m_homoclinic_orbit_coordsys;

    const ScalarType m_gain_factor { 75e-11 };
};

}

TEST(Pcr3bp_proof, homoclinic_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_homoclinic_coverings();
}

TEST(Pcr3bp_proof, periodic_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_periodic_coverings();
}

TEST(Pcr3bp_proof, jump_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_jump_coverings();
}

TEST(Pcr3bp_proof, parallelogram_coverings)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_parallelogram_coverings();
}
