///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"

#include "tools/coordsys_vec.hpp"
#include "tools/solution_curve_with_collision_check.hpp"
#include "tools/change_coordsys.hpp"
#include "tools/ry_functions.hpp"

#include "covering_relations_setup.hpp"

#include "g_map.hpp"
#include "covering_relation_checker.hpp"

#include "advanced_cone_conditions.hpp"

namespace Ursa
{

template<typename MapT>
class CoveringRelationsTest
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = Carina::LocalCoordinateSystem<MapT>;

    CoveringRelationsTest( const CoveringRelationsSetup& setup)
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

            const Carina::LocalCoordinateSystem<MapT> coordsys_src = m_homoclinic_orbit_coordsys.at(src_idx);
            const Carina::LocalCoordinateSystem<MapT> coordsys_dst = m_homoclinic_orbit_coordsys.at(dst_idx);

            const ScalarType time_span = check_covering_relation_forward(coordsys_src, coordsys_dst);
            check_collision_avoidance(coordsys_src, coordsys_dst, time_span);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check covering relations along periodic orbit
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void check_periodic_coverings()
    {
        {
            std::cout << "periodic orbit covering 0 => 1\n";

            const Carina::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(0);
            const Carina::LocalCoordinateSystem<MapT> coordsys_dst = m_periodic_orbit_coordsys.at(1);
            check_covering_relation_forward(coordsys_src, coordsys_dst);
        }

        {
            std::cout << "periodic orbit covering 1 => 2\n";

            const Carina::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(1);
            const Carina::LocalCoordinateSystem<MapT> coordsys_dst = m_periodic_orbit_coordsys.at(2);
            const ScalarType time_span = check_covering_relation_forward(coordsys_src, coordsys_dst);
            check_collision_avoidance(coordsys_src, coordsys_dst, time_span);
        }

        {
            std::cout << "periodic orbit covering 2 <= 3\n";

            const Carina::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(2);
            const Carina::LocalCoordinateSystem<MapT> coordsys_dst = m_periodic_orbit_coordsys.at(3);
            const ScalarType time_span = check_covering_relation_backward(coordsys_src, coordsys_dst);
            check_collision_avoidance(coordsys_src, coordsys_dst, time_span);
            
        }

        {
            std::cout << "periodic orbit covering 3 <= 0\n";

            const Carina::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(3);
            const Carina::LocalCoordinateSystem<MapT> coordsys_dst = m_periodic_orbit_coordsys.at(0);
            check_covering_relation_backward(coordsys_src, coordsys_dst);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check covering relations between homoclinic and periodic orbits
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void check_jump_coverings()
    {
        {
            const Carina::LocalCoordinateSystem<MapT> coordsys_src = *( m_homoclinic_orbit_coordsys.rbegin() );
            const Carina::LocalCoordinateSystem<MapT> coordsys_dst = m_periodic_orbit_coordsys.at(1);

            std::cout << "last homoclinic => periodic (1) covering\n";

            const ScalarType time_span = check_covering_relation_forward(coordsys_src, coordsys_dst);
            check_collision_avoidance(coordsys_src, coordsys_dst, time_span);
        }
        {
            const Carina::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(3);
            const Carina::LocalCoordinateSystem<MapT> coordsys_dst = *( m_homoclinic_orbit_coordsys.begin() );

            std::cout << "periodic (3) <= first homoclinic covering\n";

            const ScalarType time_span = check_covering_relation_backward(coordsys_src, coordsys_dst);
            check_collision_avoidance(coordsys_src, coordsys_dst, time_span);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check collision manifold derivative
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void collision_manifold_derivative()
    {
        LocalPoincare4_Constraint<MapT> extension_to_4_P0
        {
            std::ref(m_basic_objects.m_hamiltonian_reg2),
            std::ref( m_periodic_orbit_coordsys.at(0) )
        };

        Carina::AffineMap<MapT> P_0_c { m_periodic_orbit_coordsys.at(0) };

        using Carina::Node;
        auto collision_condition_f = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            out[0] = sqr(in[2]) + sqr(in[3]) - 8 * param[0];
        };

        MapT collision_condition(collision_condition_f, 4, 1, 1);
        collision_condition.setParameter(0, m_basic_objects.m_setup.get_mu(2));

        Carina::CompositeMap<MapT,
            decltype(extension_to_4_P0)&,
            decltype(P_0_c)&,
            decltype(collision_condition)&> composite(
                std::ref(extension_to_4_P0),
                std::ref(P_0_c),
                std::ref(collision_condition)
            );

        std::cout.precision(20);

        MatrixType der(1, 2);
        composite(N * m_gain_factor, der);
        const ScalarType fx = der(1,1);
        const ScalarType fy = der(1,2);

        std::ofstream fs("collision_manifold_parameters.txt");
        if (fs)
        {
            Carina::VariablePrinter<MapT>::print("Collision manifold derivative -ddy/ddx", -fy / fx);
            Carina::VariablePrinter<MapT>::print("Collision manifold derivative -ddx/ddy", -fx / fy);
            fs.close();
        }
        else
        {
            throw std::logic_error("Failed to export collision manifold parameters!");
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check parallelogram covering relations
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void check_parallelogram_coverings()
    {
        ScalarType alpha {};
        ScalarType p {};
        collision_avoidance_check(alpha, p);

        parallelogram_covering_check(alpha, p);
    }

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check forward covering relation between given local coordinate systems
    //! @return Time interval of underlying evolved trajectory
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ScalarType check_covering_relation_forward(
        Carina::LocalCoordinateSystem<MapT> coordsys_src,
        Carina::LocalCoordinateSystem<MapT> coordsys_dst)
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

        print_var(cr.get_der());
        print_var( Carina::span_matrix(cr.get_der()) );

        print_var(cr.get_img());
        print_var(cr.get_img_left());
        print_var(cr.get_img_right());

        EXPECT_TRUE(cr.contraction_condition());
        EXPECT_TRUE(cr.left_expansion_condition());
        EXPECT_TRUE(cr.right_expansion_condition());

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
        Carina::LocalCoordinateSystem<MapT> coordsys_src,
        Carina::LocalCoordinateSystem<MapT> coordsys_dst)
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

        MapT J = RY_Functions<MapT>::J();
        Carina::CompositeMap<MapT, decltype(J)&, decltype(f)&, decltype(J)&> jfj
        {
            std::ref(J),
            std::ref(f),
            std::ref(J)
        };

        CoveringRelationCheck cr { jfj };

        const ScalarType time_span = f.get_last_evaluation_return_time();

        print_var(cr.get_der());
        print_var( Carina::span_matrix(cr.get_der()) );

        print_var(cr.get_img());
        print_var(cr.get_img_left());
        print_var(cr.get_img_right());

        EXPECT_TRUE(cr.contraction_condition());
        EXPECT_TRUE(cr.left_expansion_condition());
        EXPECT_TRUE(cr.right_expansion_condition());

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

    void check_collision_avoidance(
        Carina::LocalCoordinateSystem<MapT> coordsys_src,
        Carina::LocalCoordinateSystem<MapT> coordsys_dst,
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
            using Carina::Node;
            auto collision_condition_f = [](Node, Node in[], int, Node out[], int, Node param[], int)
            {
                out[0] = in[0];
                out[1] = in[1];
                out[2] = sqr(in[2]) + sqr(in[3]) - 8 * param[0];
            };

            MapT collision_condition(collision_condition_f, 4, 3, 1);
            collision_condition.setParameter(0, m_basic_objects.m_setup.get_mu(2));

            Carina::MaxNorm<MapT> norm {};

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

                SolutionCurveWithCollisionCheck<MapT> solution_curve {};
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

                SolutionCurveWithCollisionCheck<MapT> solution_curve {};
                f_neg(N, time_span, solution_curve);

                bool const solution_curve_condition = solution_curve.is_condition_never_satisfied( collision_condition );
                EXPECT_TRUE(solution_curve_condition);
            }
        }
    }

    void collision_avoidance_check(ScalarType& alpha, ScalarType& p)
    {
        const VectorType arg = N;

        std::list<MatrixType> der_list {};

        MatrixType dj = {{ 0, 1 }, { 1, 0 } };

        for (int i = 0; i < 4; ++i)
        {
            const int first = i;
            const int second = (i+1) % 4;
            
            {
                G_Map<MapT> poincare
                {
                    m_basic_objects.m_vf_reg_pos2,
                    m_basic_objects.m_hamiltonian_reg2,
                    m_basic_objects.m_order,
                    m_periodic_orbit_coordsys.at(first),
                    m_periodic_orbit_coordsys.at(second),
                    m_gain_factor
                };

                MatrixType der(2,2);
                poincare(arg, der);
                print_var(der);

                der_list.emplace_back(der);
            }

            {
                G_Map<MapT> poincare
                {
                    m_basic_objects.m_vf_reg_neg2,
                    m_basic_objects.m_hamiltonian_reg2,
                    m_basic_objects.m_order,
                    m_periodic_orbit_coordsys.at(second),
                    m_periodic_orbit_coordsys.at(first),
                    m_gain_factor
                };

                MatrixType der(2,2);
                poincare(arg, der);
                der = dj * der * dj;
                print_var(der);

                der_list.emplace_back(der);
            }
        }

        auto it = der_list.begin();

        MatrixType der_union = *it;
        for (++it; it != der_list.end(); ++it)
        {
            der_union = capd::vectalg::intervalHull(*it, der_union);
        }

        print_var(der_union);
        AdvancedConeConditions<MapT> advanced_cone_conditions { der_union };

        alpha = advanced_cone_conditions.alpha;
        p = advanced_cone_conditions.p;
    }

    void parallelogram_covering_check(ScalarType alpha, ScalarType p)
    {
        const ScalarType b0 = 1.0 - 1e-7;
        const ScalarType a0 = 0.59;
        const ScalarType w0 = b0 * alpha;

        MapT R_inverse = RY_Functions<MapT>::R_Inverse(w0, a0, b0);
        MapT Y_inverse = RY_Functions<MapT>::Y_Inverse(p);
        MapT J = RY_Functions<MapT>::J();
        MapT J2 = RY_Functions<MapT>::J2();

        {
            const Carina::LocalCoordinateSystem<MapT> coordsys_src = *( m_homoclinic_orbit_coordsys.rbegin() );
            const Carina::LocalCoordinateSystem<MapT> coordsys_dst = m_periodic_orbit_coordsys.at(1);

            G_Map<MapT> poincare
            {
                m_basic_objects.m_vf_reg_pos2,
                m_basic_objects.m_hamiltonian_reg2,
                m_basic_objects.m_order,
                coordsys_src,
                coordsys_dst,
                m_gain_factor
            };

            Carina::CompositeMap<MapT,
                decltype(J2)&,
                decltype(poincare)&,
                decltype(J2)&,
                decltype(Y_inverse)&,
                decltype(R_inverse)&> composite
            {
                std::ref(J2),
                std::ref(poincare),
                std::ref(J2),
                std::ref(Y_inverse),
                std::ref(R_inverse)
            };
            
            CoveringRelationCheck cr { composite };

            print_var(cr.get_der());
            print_var( Carina::span_matrix(cr.get_der()) );

            print_var(cr.get_img());
            print_var(cr.get_img_left());
            print_var(cr.get_img_right());

            EXPECT_TRUE(cr.contraction_condition());
            EXPECT_TRUE(cr.left_expansion_condition());
            EXPECT_TRUE(cr.right_expansion_condition());
        }

        {
            const Carina::LocalCoordinateSystem<MapT> coordsys_src = m_periodic_orbit_coordsys.at(3);
            const Carina::LocalCoordinateSystem<MapT> coordsys_dst = *( m_homoclinic_orbit_coordsys.begin() );

            G_Map<MapT> poincare
            {
                m_basic_objects.m_vf_reg_neg2,
                m_basic_objects.m_hamiltonian_reg2,
                m_basic_objects.m_order,
                coordsys_dst,
                coordsys_src,
                m_gain_factor
            };

            Carina::CompositeMap<MapT,
                decltype(J)&,
                decltype(J2)&,
                decltype(poincare)&,
                decltype(J2)&,
                decltype(J)&,
                decltype(Y_inverse)&,
                decltype(R_inverse)&> composite
            {
                std::ref(J),
                std::ref(J2),
                std::ref(poincare),
                std::ref(J2),
                std::ref(J),
                std::ref(Y_inverse),
                std::ref(R_inverse)
            };
            
            CoveringRelationCheck cr { composite };

            print_var(cr.get_der());
            print_var( Carina::span_matrix(cr.get_der()) );

            print_var(cr.get_img());
            print_var(cr.get_img_left());
            print_var(cr.get_img_right());

            EXPECT_TRUE(cr.contraction_condition());
            EXPECT_TRUE(cr.left_expansion_condition());
            EXPECT_TRUE(cr.right_expansion_condition());
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
    using namespace Ursa;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_homoclinic_coverings();
}

TEST(Pcr3bp_proof, periodic_coverings)
{
    using namespace Ursa;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_periodic_coverings();
}

TEST(Pcr3bp_proof, jump_coverings)
{
    using namespace Ursa;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_jump_coverings();
}

TEST(Pcr3bp_proof, parallelogram_coverings)
{
    using namespace Ursa;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};
    CoveringRelationsTest<IMap> test { setup };
    test.check_parallelogram_coverings();
}

