///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"

#include "tools/affine_poincare_map.hpp"
#include "tools/coordsys4_alignment.hpp"
#include "tools/variable_printer.hpp"

#include "pcr3bp_reg_basic_objects.hpp"
#include "pcr3bp_reg2_initial_coordsys_generator.hpp"

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Generate local coordinate systems for homoclinic orbit
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class HomoclinicOrbitCoordsysGenerator
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = CapdUtils::LocalCoordinateSystem<MapT>;

    using AffinePoincareMap = CapdUtils::AffinePoincareMap<MapT>;

    static_assert(std::is_same<MapT, RMap>::value);

    HomoclinicOrbitCoordsysGenerator(
        const std::vector<Coordsys>& periodic_orbit_coordsys,
        const std::vector<VectorType>& homoclinic_orbit_orgins,
        ScalarType total_expansion_factor_pos,
        ScalarType total_expansion_factor_neg)
            : m_periodic_orbit_coordsys(periodic_orbit_coordsys)
            , m_homoclinic_orbit_origins(homoclinic_orbit_orgins)
    {
        const std::list<Coordsys> homoclinic_orbit_coordsys_initial = build_homoclinic_orbit_coordsys_initial();
        m_homoclinic_orbit_coordsys = build_homoclinic_orbit_coordsys(homoclinic_orbit_coordsys_initial, total_expansion_factor_pos, total_expansion_factor_neg);
    }

    const std::vector<Coordsys>& get_coordsys_container() const noexcept
    {
        return m_homoclinic_orbit_coordsys;
    }

private:
    std::list<Coordsys> build_homoclinic_orbit_coordsys_initial()
    {
        std::list<Coordsys> ret {};
        for (const auto& v : m_homoclinic_orbit_origins)
        {
            const Coordsys cs = Pcr3bp::Reg2_InitialCoordsysGenerator<MapT>::gen(m_basic_objects, v);
            ret.emplace_back(cs);
        }

        // for first coordinate system take aligned coordsys from fixed-point (w0) coordsys

        *(ret.begin()) = Coordsys(
            ret.begin()->get_origin(),
            m_periodic_orbit_coordsys.at(0).get_directions_matrix()
        );

        return ret;
    };

    std::vector<Coordsys> build_homoclinic_orbit_coordsys(
        const std::list<Coordsys>& homoclinic_orbit_coordsys_initial,
        ScalarType total_expansion_factor_pos,
        ScalarType total_expansion_factor_neg)
    {
        std::list<AffinePoincareMap> poincare_pos_list {};
        std::list<AffinePoincareMap> poincare_neg_list {};
        {
            auto it = homoclinic_orbit_coordsys_initial.begin();
            for (auto jt = std::next(it, 1); jt != homoclinic_orbit_coordsys_initial.end(); ++it, ++jt )
            {
                poincare_pos_list.emplace_back(
                    std::ref(m_basic_objects.m_vf_reg_pos2),
                    m_basic_objects.m_order,
                    *it,
                    *jt
                );

                poincare_neg_list.emplace_back(
                    std::ref(m_basic_objects.m_vf_reg_neg2),
                    m_basic_objects.m_order,
                    *jt,
                    *it
                );
            }
            poincare_neg_list.reverse();
        }

        const ScalarType expansion_factor_pos = std::pow( total_expansion_factor_pos, 0.5 / poincare_pos_list.size() );
        const ScalarType expansion_factor_neg = std::pow( total_expansion_factor_neg, 0.5 / poincare_neg_list.size() );

        CapdUtils::VariablePrinter<MapT>::print(
            "homoclinic_orbit_average_expansion_factor_pos.txt",
            "Average expansion factor along homoclinic orbit (positive direction)",
            expansion_factor_pos);

        CapdUtils::VariablePrinter<MapT>::print(
            "homoclinic_orbit_average_expansion_factor_neg.txt",
            "Average expansion factor along homoclinic orbit (negative direction)",
            expansion_factor_neg);

        const std::list<VectorType> unstable_dirs_pos = get_unstable_dirs(poincare_pos_list, VectorType{ 1.0, 0.0, 0.0, 0.0 }, expansion_factor_pos);

        const Coordsys& coordsysK = *(homoclinic_orbit_coordsys_initial.rbegin());
        const VectorType unstable_dir_pos_wK_local = *(unstable_dirs_pos.rbegin());
        const VectorType unstable_dir_pos_wK = coordsysK.get_directions_matrix() * unstable_dir_pos_wK_local;
        const VectorType stable_dir_pos_wK = AuxiliaryFunctions<MapT>::S_symmetry(unstable_dir_pos_wK);

        MatrixType coordsysK_dirs = coordsysK.get_directions_matrix();
        coordsysK_dirs.Transpose();

        const VectorType stable_dir_pos_wK_local = coordsysK_dirs * stable_dir_pos_wK;

        std::list<VectorType> unstable_dirs_neg = get_unstable_dirs(poincare_neg_list, stable_dir_pos_wK_local, expansion_factor_neg);
        unstable_dirs_neg.reverse();
        
        std::vector<Coordsys> ret {};
        ret.reserve(30);

        auto it = homoclinic_orbit_coordsys_initial.begin();
        ret.push_back( *it++ );

        size_t index = 1;

        auto it_pos = unstable_dirs_pos.begin();
        auto it_neg = std::next(unstable_dirs_neg.begin(), 1);
        for (; it != std::prev(homoclinic_orbit_coordsys_initial.end(), 1); ++it, ++it_pos, ++it_neg, ++index)
        {
            const Coordsys& cs_init = *it;
            const VectorType& p = cs_init.get_directions_matrix() * (*it_pos);
            const VectorType& n = cs_init.get_directions_matrix() * (*it_neg);

            Coordsys cs = Coordsys4_Alignment<MapT>::replace_unstable_dirs( cs_init, p, n );
            ret.push_back( cs );
        }

        {
            const Coordsys& cs_init = *it;
            const VectorType& p = cs_init.get_directions_matrix() * (*it_pos);
            const Coordsys cs = Coordsys4_Alignment<MapT>::replace_unstable_dirs_and_make_S_backsymmetric( cs_init, p );
            ret.push_back( cs );
        }

        return ret;
    }

    std::list<VectorType> get_unstable_dirs(std::list<AffinePoincareMap>& poincare_map_list, VectorType v, ScalarType expansion_factor)
    {
        std::list<VectorType> ret {};
        CapdUtils::MaxNorm<MapT> norm {};

        for (AffinePoincareMap& poincare_map : poincare_map_list)
        {
            MatrixType der(4,4);
            const VectorType val = poincare_map( VectorType(4), der );

            assert_with_exception( norm(val) < 9.7e-12 );

            auto v_init = v.euclNorm();
            v = der * v;
            v /= expansion_factor;

            ret.emplace_back( v );
        }

        return ret;
    }

    const std::vector<Coordsys>& m_periodic_orbit_coordsys;
    const std::vector<VectorType>& m_homoclinic_orbit_origins;

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    std::vector<Coordsys> m_homoclinic_orbit_coordsys {};
};

}
