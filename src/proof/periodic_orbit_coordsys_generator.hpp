///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "pcr3bp_reg_basic_objects.hpp"
#include "pcr3bp_reg2_initial_coordsys_generator.hpp"

#include <carina/poincare_wrapper.hpp>
#include <carina/timemap_wrapper.hpp>
#include <carina/composite_map.hpp>
#include <carina/extension_map.hpp>
#include <carina/projection_map.hpp>
#include <carina/affine_map.hpp>
#include <carina/constrained_function.hpp>
#include <carina/parallel_shooting/parallel_shooting_init.hpp>

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

    using Coordsys = Carina::LocalCoordinateSystem<MapT>;

    static_assert(std::is_same<MapT, RMap>::value);

    PeriodicOrbitCoordsysGenerator()
    {
        Carina::VariablePrinter<MapT>::print(
            "periodic_orbit_total_expansion_factor_pos.txt",
            "Total expansion factor along periodic orbit (positive direction)",
            m_unstable_dir_gen.get_expansion_pos_factor());

        Carina::VariablePrinter<MapT>::print(
            "periodic_orbit_total_expansion_factor_neg.txt",
            "Total expansion factor along periodic orbit (negative direction)",
            m_unstable_dir_gen.get_expansion_pos_factor());

        Carina::VariablePrinter<MapT>::print(
            "periodic_orbit_average_expansion_factor_pos.txt",
            "Average expansion factor along periodic orbit (positive direction)",
            m_expansion_factor_pos);

        Carina::VariablePrinter<MapT>::print(
            "periodic_orbit_average_expansion_factor_neg.txt",
            "Average expansion factor along periodic orbit (negative direction)",
            m_expansion_factor_neg);

        m_local_poincare_pos.at(0)(VectorType(2));

        Carina::VariablePrinter<MapT>::print(
            "periodic_orbit_g_0_1_approx_return_time.txt",
            "Approximate value of return time on underlying Poincare map of g_01 map",
            m_local_poincare_pos.at(0).get_last_evaluation_return_time() );

        m_local_poincare_pos.at(1)(VectorType(2));

        Carina::VariablePrinter<MapT>::print(
            "periodic_orbit_g_1_2_approx_return_time.txt",
            "Approximate value of return time on underlying Poincare map of g_12 map",
            m_local_poincare_pos.at(1).get_last_evaluation_return_time() );
    }

    const std::vector<Coordsys>& get_coordsys_container() const
    {
        return m_local_coord;
    }

private:
    class UnstableDirsList : public std::list<VectorType>
    {
    public:
        UnstableDirsList(const std::list<VectorType>& arg)
        {
            assert_with_exception(arg.size() == 3);
            std::list<VectorType>::operator= (arg);
        }
    };

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};
    RegLyapunovCollisionOrbitParameters<MapT> m_parameters {};

    std::array<Coordsys, 4> m_initial_coordsys
    {
        Pcr3bpReg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_initial_point()),
        Pcr3bpReg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_intermediate_point()),
        Pcr3bpReg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_image_point()),
        Pcr3bpReg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_intermediate_point_neg())
    };

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

    Carina::CompositeMap<MapT, LocalPoincare4<MapT>&, LocalPoincare4<MapT>&> m_local_poincare_pos_dual
    {
        std::ref(m_local_poincare_pos.at(0)),
        std::ref(m_local_poincare_pos.at(1))
    };

    Carina::CompositeMap<MapT, LocalPoincare4<MapT>&, LocalPoincare4<MapT>&> m_local_poincare_neg_dual
    {
        std::ref(m_local_poincare_neg.at(0)),
        std::ref(m_local_poincare_neg.at(1))
    };

    Carina::UnstableDirectionsGenerator<MapT,
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

    Carina::DirectionShifting<MapT,
        LocalPoincare4<MapT>&,
        LocalPoincare4<MapT>&,
        LocalPoincare4<MapT>&> m_direction_shifting_pos
    {
        m_expansion_factor_pos,
        std::ref(m_local_poincare_pos.at(0)),
        std::ref(m_local_poincare_pos.at(1)),
        std::ref(m_local_poincare_pos.at(2))
    };

    Carina::DirectionShifting<MapT,
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
    };

    const std::vector<Coordsys> m_local_coord
    {
        Coordsys4_Alignment<MapT>::align_with_s_symmetry(
            m_initial_coordsys.at(0),
            m_unstable_dir_gen.get_unstable_pos()),
        Coordsys4_Alignment<MapT>::align(
            m_initial_coordsys.at(1),
            *std::next(m_unstable_pos_dirs_2d.begin(), 0),
            *std::next(m_unstable_neg_dirs_2d.rbegin(), 0)),
        Coordsys4_Alignment<MapT>::align(
            m_initial_coordsys.at(2),
            *std::next(m_unstable_pos_dirs_2d.begin(), 1),
            *std::next(m_unstable_neg_dirs_2d.rbegin(), 1)),
        Coordsys4_Alignment<MapT>::align(
            m_initial_coordsys.at(3),
            *std::next(m_unstable_pos_dirs_2d.begin(), 2),
            *std::next(m_unstable_neg_dirs_2d.rbegin(), 2))
    };
};

}
