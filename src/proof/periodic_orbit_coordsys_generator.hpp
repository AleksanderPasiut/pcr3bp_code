///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "pcr3bp_reg_basic_objects.hpp"
#include "pcr3bp_reg2_initial_coordsys_generator.hpp"

#include <capd_utils/poincare_wrapper.hpp>
#include <capd_utils/timemap_wrapper.hpp>
#include <capd_utils/composite_map.hpp>
#include <capd_utils/extension_map.hpp>
#include <capd_utils/projection_map.hpp>
#include <capd_utils/affine_map.hpp>
#include <capd_utils/constrained_function.hpp>
#include <capd_utils/parallel_shooting/parallel_shooting_init.hpp>

#include "tools/local_poincare4.hpp"
#include "tools/coordsys4_alignment.hpp"
#include "tools/unstable_directions_generator.hpp"
#include "tools/direction_shifting.hpp"
#include "tools/auxiliary_functions.hpp"
#include "tools/variable_printer.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Compute coordinate systems in fixed point and in 3 other points that are approximately located on periodic orbit.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class PeriodicOrbitCoordsysGenerator
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = CapdUtils::LocalCoordinateSystem<MapT>;

    static_assert(std::is_same<MapT, RMap>::value);

    PeriodicOrbitCoordsysGenerator()
    {
        std::cout.precision(15);
        // {
        //     print_var( m_local_coord[0].get_directions_matrix() );
        //     print_var( m_local_coord[1].get_directions_matrix() );
        //     print_var( m_local_coord[2].get_directions_matrix() );
        //     print_var( m_local_coord[3].get_directions_matrix() );
        // }
        
        {
            CapdUtils::AffinePoincareMap poincare_1_pos
            {
                m_basic_objects.m_vf_reg_pos2,
                m_basic_objects.m_order,
                m_initial_coordsys.at(0),
                m_initial_coordsys.at(1)
            };

            CapdUtils::AffinePoincareMap poincare_2_pos
            {
                m_basic_objects.m_vf_reg_pos2,
                m_basic_objects.m_order,
                m_initial_coordsys.at(1),
                m_initial_coordsys.at(2)
            };

            CapdUtils::AffinePoincareMap poincare_3_pos
            {
                m_basic_objects.m_vf_reg_pos2,
                m_basic_objects.m_order,
                m_initial_coordsys.at(2),
                m_initial_coordsys.at(3)
            };
            
            CapdUtils::AffinePoincareMap poincare_0_pos
            {
                m_basic_objects.m_vf_reg_pos2,
                m_basic_objects.m_order,
                m_initial_coordsys.at(3),
                m_initial_coordsys.at(0)
            };
            
            CapdUtils::AffinePoincareMap poincare_1_neg
            {
                m_basic_objects.m_vf_reg_neg2,
                m_basic_objects.m_order,
                m_initial_coordsys.at(2),
                m_initial_coordsys.at(1)
            };

            CapdUtils::AffinePoincareMap poincare_0_neg
            {
                m_basic_objects.m_vf_reg_neg2,
                m_basic_objects.m_order,
                m_initial_coordsys.at(1),
                m_initial_coordsys.at(0)
            };
            
            CapdUtils::CompositeMap<MapT,
                CapdUtils::AffinePoincareMap<MapT>&,
                CapdUtils::AffinePoincareMap<MapT>&,
                CapdUtils::AffinePoincareMap<MapT>&,
                CapdUtils::AffinePoincareMap<MapT>&> poincare_total
            {
                std::ref(poincare_1_pos),
                std::ref(poincare_2_pos),
                std::ref(poincare_3_pos),
                std::ref(poincare_0_pos)
            };

            const VectorType zero4 = VectorType{ 0.0, 0.0, 0.0, 0.0 };

            MatrixType der {};
            auto x1 = poincare_total( zero4, der );
            print_var(x1);

            print_var( der );

            const VectorType unstable_dir_w0_local = CapdUtils::PowerIteration<MapT>::evaluate( der, VectorType{ 1.0, 0.0, 0.0, 0.0 }, 100 );
            print_var( unstable_dir_w0_local );
            print_var( der * unstable_dir_w0_local );

            print_var( unstable_dir_w0_local.euclNorm() );
            print_var( (der * unstable_dir_w0_local).euclNorm() );

            const ScalarType expansion_factor = std::pow( (der * unstable_dir_w0_local).euclNorm(), 0.25 );
            print_var(expansion_factor);

            const VectorType unstable_dir_w0 = m_initial_coordsys.at(0).get_directions_matrix() * unstable_dir_w0_local;
            print_var( unstable_dir_w0 );

            // print_var( m_unstable_dir_gen.get_expansion_pos_factor() ); 
            const VectorType stable_dir_w0 = AuxiliaryFunctions<MapT>::S_symmetry(unstable_dir_w0);
            print_var( stable_dir_w0 );
            

            MatrixType der1 {};
            print_var( poincare_1_pos(zero4, der1) );
            const VectorType unstable_dir_w1_local = (der1 * unstable_dir_w0_local) / expansion_factor;
            const VectorType unstable_dir_w1 = m_initial_coordsys.at(1).get_directions_matrix() * unstable_dir_w1_local;
            print_var( unstable_dir_w1 );

            MatrixType der2 {};
            print_var( poincare_2_pos(zero4, der2) );
            const VectorType unstable_dir_w2_local = (der2 * unstable_dir_w1_local) / expansion_factor;
            const VectorType unstable_dir_w2 = m_initial_coordsys.at(2).get_directions_matrix() * unstable_dir_w2_local;
            print_var( unstable_dir_w2 );

            const VectorType stable_dir_w2 = AuxiliaryFunctions<MapT>::S_symmetry(unstable_dir_w2);
            print_var( stable_dir_w2 );

            MatrixType w2_initial_dirs = m_initial_coordsys.at(2).get_directions_matrix();
            w2_initial_dirs.Transpose();

            const VectorType stable_dir_w2_local = w2_initial_dirs * stable_dir_w2;
            print_var( stable_dir_w2_local );

            MatrixType der1_neg {};
            print_var( poincare_1_neg(zero4, der1_neg) );
            const VectorType stable_dir_w1_local = (der1_neg * stable_dir_w2_local) / expansion_factor;
            const VectorType stable_dir_w1 = m_initial_coordsys.at(1).get_directions_matrix() * stable_dir_w1_local;
            print_var( stable_dir_w1 );

            m_local_coord = std::vector<Coordsys>
            {
                Coordsys4_Alignment<MapT>::replace_unstable_dirs(
                    m_initial_coordsys.at(0),
                    unstable_dir_w0,
                    stable_dir_w0),
                Coordsys4_Alignment<MapT>::replace_unstable_dirs(
                    m_initial_coordsys.at(1),
                    unstable_dir_w1,
                    stable_dir_w1),
                Coordsys4_Alignment<MapT>::replace_unstable_dirs(
                    m_initial_coordsys.at(2),
                    unstable_dir_w2,
                    stable_dir_w2),
                Coordsys4_Alignment<MapT>::replace_unstable_dirs(
                    m_initial_coordsys.at(3),
                    AuxiliaryFunctions<MapT>::S_symmetry(stable_dir_w1),
                    AuxiliaryFunctions<MapT>::S_symmetry(unstable_dir_w1))
            };




            // VectorType x0 = m_local_coord[0].get_origin();
            // print_var( x0 );

            // MatrixType der_pos(4, 4);
            // timemap_pos( x0, der_pos );

            // MatrixType der_neg(4, 4);
            // timemap_neg( x0, der_neg );

            // const MatrixType der = CapdUtils::gaussInverseMatrix<MapT>(der_neg) * der_pos;

            // const VectorType v3 = CapdUtils::Extract<MapT>::get_vvector(dirs, 3);
            // const VectorType v4 = CapdUtils::Extract<MapT>::get_vvector(dirs, 4);

            // const VectorType unstable_dir = CapdUtils::PowerIteration<MapT>::evaluate( der, VectorType{ 1.0, 0.0, 0.0, 0.0 }, 100 );
            // print_var(unstable_dir);

            // const VectorType unstable_dir_projected = unstable_dir - v3 * capd::vectalg::scalarProduct(v3, unstable_dir) - v4 * capd::vectalg::scalarProduct(v4, unstable_dir);
            // print_var(unstable_dir_projected);
        }


        

        // CapdUtils::VariablePrinter<MapT>::print(
        //     "periodic_orbit_total_expansion_factor_pos.txt",
        //     "Total expansion factor along periodic orbit (positive direction)",
        //     m_unstable_dir_gen.get_expansion_pos_factor());

        // CapdUtils::VariablePrinter<MapT>::print(
        //     "periodic_orbit_total_expansion_factor_neg.txt",
        //     "Total expansion factor along periodic orbit (negative direction)",
        //     m_unstable_dir_gen.get_expansion_pos_factor());

        // CapdUtils::VariablePrinter<MapT>::print(
        //     "periodic_orbit_average_expansion_factor_pos.txt",
        //     "Average expansion factor along periodic orbit (positive direction)",
        //     m_expansion_factor_pos);

        // CapdUtils::VariablePrinter<MapT>::print(
        //     "periodic_orbit_average_expansion_factor_neg.txt",
        //     "Average expansion factor along periodic orbit (negative direction)",
        //     m_expansion_factor_neg);

        // m_local_poincare_pos.at(0)(VectorType(2));

        // CapdUtils::VariablePrinter<MapT>::print(
        //     "periodic_orbit_g_0_1_approx_return_time.txt",
        //     "Approximate value of return time on underlying Poincare map of g_01 map",
        //     m_local_poincare_pos.at(0).get_last_evaluation_return_time() );

        // m_local_poincare_pos.at(1)(VectorType(2));

        // CapdUtils::VariablePrinter<MapT>::print(
        //     "periodic_orbit_g_1_2_approx_return_time.txt",
        //     "Approximate value of return time on underlying Poincare map of g_12 map",
        //     m_local_poincare_pos.at(1).get_last_evaluation_return_time() );
    }

    const std::vector<Coordsys>& get_coordsys_container() const
    {
        return m_local_coord;
    }

private:
    /*class UnstableDirsList : public std::list<VectorType>
    {
    public:
        UnstableDirsList(const std::list<VectorType>& arg)
        {
            assert_with_exception(arg.size() == 3);
            std::list<VectorType>::operator= (arg);
        }
    };*/

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};
    RegLyapunovCollisionOrbitParameters<MapT> m_parameters {};

    std::array<Coordsys, 4> m_initial_coordsys
    {
        Pcr3bp::Reg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_initial_point()),
        Pcr3bp::Reg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_intermediate_point()),
        Pcr3bp::Reg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_image_point()),
        Pcr3bp::Reg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_intermediate_point_neg())
    };

    /*
    std::array<LocalPoincare4<MapT>, 3> m_local_poincare_pos
    {
        LocalPoincare4<MapT>
        {
            m_basic_objects.m_vf_reg_pos2,
            m_basic_objects.m_hamiltonian_reg2,
            m_basic_objects.m_order,
            m_initial_coordsys.at(0),
            m_initial_coordsys.at(1)
        },
        LocalPoincare4<MapT>
        {
            m_basic_objects.m_vf_reg_pos2,
            m_basic_objects.m_hamiltonian_reg2,
            m_basic_objects.m_order,
            m_initial_coordsys.at(1),
            m_initial_coordsys.at(2)
        },
        LocalPoincare4<MapT>
        {
            m_basic_objects.m_vf_reg_pos2,
            m_basic_objects.m_hamiltonian_reg2,
            m_basic_objects.m_order,
            m_initial_coordsys.at(2),
            m_initial_coordsys.at(3)
        }
    };

    std::array<LocalPoincare4<MapT>, 3> m_local_poincare_neg
    {
        LocalPoincare4<MapT>
        {
            m_basic_objects.m_vf_reg_neg2,
            m_basic_objects.m_hamiltonian_reg2,
            m_basic_objects.m_order,
            m_initial_coordsys.at(0),
            m_initial_coordsys.at(3)
        },
        LocalPoincare4<MapT>
        {
            m_basic_objects.m_vf_reg_neg2,
            m_basic_objects.m_hamiltonian_reg2,
            m_basic_objects.m_order,
            m_initial_coordsys.at(3),
            m_initial_coordsys.at(2)
        },
        LocalPoincare4<MapT>
        {
            m_basic_objects.m_vf_reg_neg2,
            m_basic_objects.m_hamiltonian_reg2,
            m_basic_objects.m_order,
            m_initial_coordsys.at(2),
            m_initial_coordsys.at(1)
        }
    };

    CapdUtils::CompositeMap<MapT, LocalPoincare4<MapT>&, LocalPoincare4<MapT>&> m_local_poincare_pos_dual
    {
        std::ref(m_local_poincare_pos.at(0)),
        std::ref(m_local_poincare_pos.at(1))
    };

    CapdUtils::CompositeMap<MapT, LocalPoincare4<MapT>&, LocalPoincare4<MapT>&> m_local_poincare_neg_dual
    {
        std::ref(m_local_poincare_neg.at(0)),
        std::ref(m_local_poincare_neg.at(1))
    };

    CapdUtils::UnstableDirectionsGenerator<MapT,
        decltype(m_local_poincare_pos_dual)&,
        decltype(m_local_poincare_neg_dual)&> m_unstable_dir_gen
    {
        std::ref(m_local_poincare_pos_dual),
        std::ref(m_local_poincare_neg_dual),
        VectorType{ 1.0, 0.0 },
        50
    };

    ScalarType const m_expansion_factor_pos
    {
        std::pow(static_cast<double>(m_unstable_dir_gen.get_expansion_pos_factor()), 0.25)
    };

    ScalarType const m_expansion_factor_neg
    {
        std::pow(static_cast<double>(m_unstable_dir_gen.get_expansion_neg_factor()), 0.25)
    };

    CapdUtils::DirectionShifting<MapT,
        LocalPoincare4<MapT>&,
        LocalPoincare4<MapT>&,
        LocalPoincare4<MapT>&> m_direction_shifting_pos
    {
        m_expansion_factor_pos,
        std::ref(m_local_poincare_pos.at(0)),
        std::ref(m_local_poincare_pos.at(1)),
        std::ref(m_local_poincare_pos.at(2))
    };

    CapdUtils::DirectionShifting<MapT,
        LocalPoincare4<MapT>&,
        LocalPoincare4<MapT>&,
        LocalPoincare4<MapT>&> m_direction_shifting_neg
    {
        m_expansion_factor_neg,
        std::ref(m_local_poincare_neg.at(0)),
        std::ref(m_local_poincare_neg.at(1)),
        std::ref(m_local_poincare_neg.at(2))
    };

    const UnstableDirsList m_unstable_pos_dirs_2d 
    {
        m_direction_shifting_pos.eval( m_unstable_dir_gen.get_unstable_pos(), { VectorType(2), VectorType(2), VectorType(2) } )
    };

    const UnstableDirsList m_unstable_neg_dirs_2d
    {
        m_direction_shifting_neg.eval( m_unstable_dir_gen.get_unstable_neg(), { VectorType(2), VectorType(2), VectorType(2) } )
    };*/

    std::vector<Coordsys> m_local_coord {};
    // {
    //     Coordsys4_Alignment<MapT>::align_with_s_symmetry(
    //         m_initial_coordsys.at(0),
    //         m_unstable_dir_gen.get_unstable_pos()),
    //     Coordsys4_Alignment<MapT>::align(
    //         m_initial_coordsys.at(1),
    //         *std::next(m_unstable_pos_dirs_2d.begin(), 0),
    //         *std::next(m_unstable_neg_dirs_2d.rbegin(), 0)),
    //     Coordsys4_Alignment<MapT>::align_with_s_symmetry(
    //         m_initial_coordsys.at(2),
    //         *std::next(m_unstable_pos_dirs_2d.begin(), 1)),
    //     Coordsys4_Alignment<MapT>::align(
    //         m_initial_coordsys.at(3),
    //         *std::next(m_unstable_pos_dirs_2d.begin(), 2),
    //         *std::next(m_unstable_neg_dirs_2d.rbegin(), 2))
    // };
};

}
