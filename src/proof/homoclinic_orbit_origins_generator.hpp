///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "pcr3bp_reg_basic_objects.hpp"

#include <carina/poincare_wrapper.hpp>
#include <carina/timemap_wrapper.hpp>

#include "homoclinic_orbit_origins_initial.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Generate homoclinic orbit origins
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class HomoclinicOrbitOriginsGenerator
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    HomoclinicOrbitOriginsGenerator(const HomoclinicOrbitOriginsInitial<MapT>& homoclinic_orbit_origins_initial)
        : m_homoclinic_orbit_origins_initial(homoclinic_orbit_origins_initial)
    {
        const std::vector<VectorType>& initial_origins = homoclinic_orbit_origins_initial.get_points();

        const double t_inter_1 = 0.895631;
        const double t_inter_2 = 0.0133115509;
        const double t_inter_3 = 1.041086989;

        auto move_pt = [this]( VectorType pt, ScalarType time) -> VectorType
        {
            if (time > 0.0)
            {
                m_timemap_pos.set_time(+time);
                return m_timemap_pos( pt );
            }
            else
            {
                m_timemap_neg.set_time(-time);
                return m_timemap_neg( pt );
            }
        };

        m_points.reserve(30);
        m_points.push_back( initial_origins.at(0) ); // ---------------------- c4
        m_points.push_back( move_pt( initial_origins.at(1), -t_inter_2 ) );
        m_points.push_back( initial_origins.at(1) ); // ---------------------- c6
        m_points.push_back( move_pt( initial_origins.at(1), +t_inter_2 ) );
        m_points.push_back( initial_origins.at(2) ); // ---------------------- c8
        m_points.push_back( move_pt( initial_origins.at(3), -t_inter_2 ) );
        m_points.push_back( initial_origins.at(3) ); // ---------------------- c10
        m_points.push_back( move_pt( initial_origins.at(3), +t_inter_2 ) );
        m_points.push_back( initial_origins.at(4) ); // ---------------------- c12
        m_points.push_back( move_pt( initial_origins.at(5), -t_inter_2 ) );
        m_points.push_back( initial_origins.at(5) ); // ---------------------- c14
        m_points.push_back( move_pt( initial_origins.at(5), +t_inter_2 ) );
        m_points.push_back( move_pt( initial_origins.at(5), +t_inter_2 + t_inter_3/2 ) );
        m_points.push_back( initial_origins.at(6) ); // ---------------------- c17
        m_points.push_back( initial_origins.at(7) ); // ---------------------- c18
        m_points.push_back( initial_origins.at(8) ); // ---------------------- c19
        m_points.push_back( move_pt( initial_origins.at(9), -t_inter_2 - t_inter_3/2 ) );
        m_points.push_back( move_pt( initial_origins.at(9), -t_inter_2 ) );
        m_points.push_back( initial_origins.at(9) ); // ---------------------- c22
        m_points.push_back( move_pt( initial_origins.at(9), +t_inter_2 ) );
        m_points.push_back( initial_origins.at(10) ); // ---------------------- c24
        m_points.push_back( move_pt( initial_origins.at(11), -t_inter_2 ) );
        m_points.push_back( initial_origins.at(11) ); // ---------------------- c26
        m_points.push_back( move_pt( initial_origins.at(11), +t_inter_2 ) );
        m_points.push_back( initial_origins.at(12) ); // ---------------------- c28
        m_points.push_back( move_pt( initial_origins.at(13), -t_inter_2 ) );
        m_points.push_back( initial_origins.at(13) ); // ---------------------- c30
        m_points.push_back( move_pt( initial_origins.at(13), +t_inter_2 ) );
        m_points.push_back( initial_origins.at(14) ); // ---------------------- c32
    }

    const std::vector<VectorType>& get_points() const noexcept
    {
        return m_points;
    }

    ScalarType get_total_expansion_factor_pos() const noexcept
    {
        return m_total_expansion_factor_pos;
    }

    ScalarType get_total_expansion_factor_neg() const noexcept
    {
        return m_total_expansion_factor_neg;
    }

private:
    const HomoclinicOrbitOriginsInitial<MapT>& m_homoclinic_orbit_origins_initial;

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    Carina::CoordinateSection<MapT> m_v_section { 4, 1, ScalarType(0.0) };
    Carina::PoincareWrapper<MapT, decltype(m_v_section)> m_poincare_pos
    {
        m_basic_objects.m_vf_reg_pos2,
        m_basic_objects.m_order,
        m_v_section
    };

    Carina::TimemapWrapper<MapT> m_timemap_pos
    {
        m_basic_objects.m_vf_reg_pos2,
        0.0,
        m_basic_objects.m_order
    };

    Carina::TimemapWrapper<MapT> m_timemap_neg
    {
        m_basic_objects.m_vf_reg_neg2,
        0.0,
        m_basic_objects.m_order
    };

    std::vector<VectorType> m_points {};

    const ScalarType m_total_expansion_factor_pos
    {
        m_homoclinic_orbit_origins_initial.get_total_expansion_factor_pos()
    };

    const ScalarType m_total_expansion_factor_neg
    {
        m_homoclinic_orbit_origins_initial.get_total_expansion_factor_neg()
    };
};

}
