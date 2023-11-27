///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"

#include "tools/local_poincare4.hpp"
#include "tools/coordsys4_alignment.hpp"
#include "tools/variable_printer.hpp"

#include "pcr3bp_reg_basic_objects.hpp"
#include "pcr3bp_reg2_initial_coordsys_generator.hpp"

namespace Ursa
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

    using Coordsys = Carina::LocalCoordinateSystem<MapT>;

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

        // for first and last coordinate system take aligned coordsys from fixed-point coordsys

        *(ret.begin()) = Coordsys(
            ret.begin()->get_origin(),
            m_periodic_orbit_coordsys.at(0).get_directions_matrix()
        );

        *(ret.rbegin()) = Coordsys(
            ret.rbegin()->get_origin(),
            m_periodic_orbit_coordsys.at(0).get_directions_matrix()
        );

        return ret;
    };

    std::vector<Coordsys> build_homoclinic_orbit_coordsys(
        const std::list<Coordsys>& homoclinic_orbit_coordsys_initial,
        ScalarType total_expansion_factor_pos,
        ScalarType total_expansion_factor_neg)
    {
        std::list<LocalPoincare4<MapT>> poincare_pos_list {};
        {
            auto it = homoclinic_orbit_coordsys_initial.begin();
            for (auto jt = std::next(it, 1); jt != homoclinic_orbit_coordsys_initial.end(); ++it, ++jt )
            {
                poincare_pos_list.emplace_back(
                    std::ref(m_basic_objects.m_vf_reg_pos2),
                    std::ref(m_basic_objects.m_hamiltonian_reg2),
                    m_basic_objects.m_order,
                    *it,
                    *jt
                );
            }
        }

        std::list<LocalPoincare4<MapT>> poincare_neg_list {};
        {
            auto it = homoclinic_orbit_coordsys_initial.rbegin();
            for (auto jt = std::next(it, 1); jt != homoclinic_orbit_coordsys_initial.rend(); ++it, ++jt )
            {
                poincare_neg_list.emplace_back(
                    std::ref(m_basic_objects.m_vf_reg_neg2),
                    std::ref(m_basic_objects.m_hamiltonian_reg2),
                    m_basic_objects.m_order,
                    *it,
                    *jt
                );
            }
        }

        const ScalarType expansion_factor_pos = std::pow( total_expansion_factor_pos, 1.0 / poincare_pos_list.size() );
        const ScalarType expansion_factor_neg = std::pow( total_expansion_factor_neg, 1.0 / poincare_neg_list.size() );

        Carina::VariablePrinter<MapT>::print(
            "homoclinic_orbit_average_expansion_factor_pos.txt",
            "Average expansion factor along homoclinic orbit (positive direction)",
            expansion_factor_pos);

        Carina::VariablePrinter<MapT>::print(
            "homoclinic_orbit_average_expansion_factor_neg.txt",
            "Average expansion factor along homoclinic orbit (negative direction)",
            expansion_factor_neg);

        const std::list<VectorType> unstable_dirs_pos_2d = get_unstable_dirs(poincare_pos_list, VectorType{ 1.0, 0.0 }, expansion_factor_pos);
        const std::list<VectorType> unstable_dirs_neg_2d = get_unstable_dirs(poincare_neg_list, VectorType{ 0.0, 1.0 }, expansion_factor_neg);
        
        std::vector<Coordsys> ret {};
        ret.reserve(30);

        auto it = homoclinic_orbit_coordsys_initial.begin();
        ret.push_back( *it++ );

        size_t index = 1;

        auto it_pos = unstable_dirs_pos_2d.begin();
        auto it_neg = std::next(unstable_dirs_neg_2d.rbegin(), 1);
        for (; it != std::prev(homoclinic_orbit_coordsys_initial.end(), 1); ++it, ++it_pos, ++it_neg, ++index)
        {
            const Coordsys& cs_init = *it;
            const VectorType& p = *it_pos;
            const VectorType& n = *it_neg;

            if (index == 14)
            {
                Coordsys cs = Coordsys4_Alignment<MapT>::align_with_s_symmetry( cs_init, p );
                ret.push_back( cs );
            }
            else
            {
                Coordsys cs = Coordsys4_Alignment<MapT>::align( cs_init, p, n );
                ret.push_back( cs );
            }
        }

        ret.push_back( *it );

        return ret;
    }

    std::list<VectorType> get_unstable_dirs(std::list<LocalPoincare4<MapT>>& poincare_map_list, VectorType v, ScalarType expansion_factor)
    {
        std::list<VectorType> ret {};
        Carina::MaxNorm<MapT> norm {};

        for (LocalPoincare4<MapT>& poincare_map : poincare_map_list)
        {
            MatrixType der(2,2);
            const VectorType val = poincare_map( VectorType(2), der );

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
