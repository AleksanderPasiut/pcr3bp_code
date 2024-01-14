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

#include "tools/affine_poincare_map.hpp"
#include "tools/coordsys4_alignment.hpp"
#include "tools/power_iteration.hpp"
#include "tools/auxiliary_functions.hpp"
#include "tools/variable_printer.hpp"

namespace Pcr3bpProof
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

            CapdUtils::MaxNorm<MapT> norm {};

            MatrixType der {};
            {
                const ScalarType epsilon = norm( poincare_total( VectorType(4), der ) );
                if (epsilon > 1.4e-12)
                {
                    std::cout << "WARNING at line " << __LINE__ << ": Result norm exceeds threshold! (epsilon = " << epsilon << ")";
                }
            }

            const VectorType unstable_dir_w0_local = CapdUtils::PowerIteration<MapT>::evaluate( der, VectorType{ 1.0, 0.0, 0.0, 0.0 }, 100 );
            const ScalarType expansion_factor = std::pow( (der * unstable_dir_w0_local).euclNorm(), 0.25 );

            const VectorType unstable_dir_w0 = m_initial_coordsys.at(0).get_directions_matrix() * unstable_dir_w0_local;
            const VectorType stable_dir_w0 = AuxiliaryFunctions<MapT>::S_symmetry(unstable_dir_w0);

            MatrixType der1 {};
            {
                const ScalarType epsilon = norm( poincare_1_pos(VectorType(4), der1) );
                if (epsilon > 2.9e-15)
                {
                    std::cout << "WARNING at line " << __LINE__ << ": Result norm exceeds threshold! (epsilon = " << epsilon << ")\n";
                }
            }

            const VectorType unstable_dir_w1_local = (der1 * unstable_dir_w0_local) / expansion_factor;
            const VectorType unstable_dir_w1 = m_initial_coordsys.at(1).get_directions_matrix() * unstable_dir_w1_local;

            MatrixType der2 {};
            {
                const ScalarType epsilon = norm( poincare_2_pos(VectorType(4), der2) );
                if (epsilon > 1.2e-16)
                {
                    std::cout << "WARNING at line " << __LINE__ << ": Result norm exceeds threshold! (epsilon = " << epsilon << ")\n";
                }
            }

            const VectorType unstable_dir_w2_local = (der2 * unstable_dir_w1_local) / expansion_factor;
            const VectorType unstable_dir_w2 = m_initial_coordsys.at(2).get_directions_matrix() * unstable_dir_w2_local;

            const VectorType stable_dir_w2 = AuxiliaryFunctions<MapT>::S_symmetry(unstable_dir_w2);

            MatrixType w2_initial_dirs = m_initial_coordsys.at(2).get_directions_matrix();
            w2_initial_dirs.Transpose();

            const VectorType stable_dir_w2_local = w2_initial_dirs * stable_dir_w2;

            MatrixType der1_neg {};
            {
                const ScalarType epsilon = norm( poincare_1_neg(VectorType(4), der1_neg) );
                if (epsilon > 7.2e-14)
                {
                    std::cout << "WARNING at line " << __LINE__ << ": Result norm exceeds threshold! (epsilon = " << epsilon << ")\n";
                }
            }

            const VectorType stable_dir_w1_local = (der1_neg * stable_dir_w2_local) / expansion_factor;
            const VectorType stable_dir_w1 = m_initial_coordsys.at(1).get_directions_matrix() * stable_dir_w1_local;

            m_local_coord.reserve(4);
            m_local_coord.push_back(
                Coordsys4_Alignment<MapT>::replace_unstable_dirs_and_make_S_backsymmetric(
                    m_initial_coordsys.at(0),
                    unstable_dir_w0));
            
            m_local_coord.push_back(
                Coordsys4_Alignment<MapT>::replace_unstable_dirs(
                    m_initial_coordsys.at(1),
                    unstable_dir_w1,
                    stable_dir_w1));

            m_local_coord.push_back(
                Coordsys4_Alignment<MapT>::replace_unstable_dirs_and_make_S_backsymmetric(
                    m_initial_coordsys.at(2),
                    unstable_dir_w2));

            m_local_coord.push_back(
                Coordsys4_Alignment<MapT>::create_S_backsymmetric(
                    m_local_coord.at(1)));
        }
    }

    const std::vector<Coordsys>& get_coordsys_container() const
    {
        return m_local_coord;
    }

private:
    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};
    RegLyapunovCollisionOrbitParameters<MapT> m_parameters {};

    std::array<Coordsys, 4> m_initial_coordsys
    {
        Pcr3bp::Reg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_initial_point()),
        Pcr3bp::Reg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_intermediate_point()),
        Pcr3bp::Reg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_image_point()),
        Pcr3bp::Reg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, m_basic_objects.m_parameters.get_intermediate_point_neg())
    };

    std::vector<Coordsys> m_local_coord {};
};

}
