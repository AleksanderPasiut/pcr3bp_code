///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/local_coordinate_system.hpp>
#include <capd_utils/gauss.hpp>
#include <capd_utils/extract.hpp>

#include "psi0_specialized.hpp"

#include <proof/pcr3bp_reg_basic_objects.hpp>
#include <proof/periodic_orbit_coordsys_generator.hpp>

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Compute coefficients d1 and d2 for psi0 map
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Psi0_Coefficients
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static Psi0_Coefficients& get()
    {
        static Psi0_Coefficients s_instance {};
        return s_instance;
    }
    
    MapT& get_internal_map_ref() noexcept
    {
        return m_internal_map;
    }

    const std::array<ScalarType, 2>& get_d_coeffs() const noexcept
    {
        return m_d;
    }

private:
    Psi0_Coefficients()
    {}

    static std::array<ScalarType, 2> compute_d_coeffs(
        MapT& internal_map,
        const CapdUtils::LocalCoordinateSystem<MapT>& dst_coordsys)
    {
        const VectorType unstable_dir = CapdUtils::Extract<MapT>::get_vvector(dst_coordsys.get_directions_matrix(), 1);

        MatrixType dd(4, 2);
        internal_map( VectorType(2), dd );

        MatrixType dd2(2,4);
        for (int i = 1; i <= 2; ++i)
        {
            for (int j = 1; j <= 4; ++j)
            {
                dd2(i,j) = dd(j,i);
            }
        }

        const MatrixType ddx = dd2 * dd;

        const VectorType unstable_dir_reg = dd2 * unstable_dir;
        const VectorType d_coeff = CapdUtils::gauss<MapT>(ddx, unstable_dir_reg);

        return std::array<ScalarType, 2>
        {
            d_coeff[0],
            d_coeff[1]
        };
    }

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    MapT m_internal_map
    {
        Psi0_specialized<MapT>::create( m_basic_objects.m_h0, m_basic_objects.m_setup )
    };

    PeriodicOrbitCoordsysGenerator<RMap> m_periodic_orbit_coordsys_generator_approx {};

    std::vector<CapdUtils::LocalCoordinateSystem<RMap>> m_periodic_orbit_coordsys_approx
    {
        m_periodic_orbit_coordsys_generator_approx.get_coordsys_container()
    };

    RegLyapunovCollisionOrbitParameters<IMap> m_parameters {};

    const CapdUtils::LocalCoordinateSystem<MapT> m_src_coordsys_4_dim
    {
        CapdUtils::LocalCoordinateSystem<IMap>(
            m_parameters.get_initial_point(),
            CapdUtils::matrix_cast<IMatrix>(m_periodic_orbit_coordsys_approx.at(0).get_directions_matrix()) )
    };

    const std::array<ScalarType, 2> m_d
    {
        compute_d_coeffs(m_internal_map, m_src_coordsys_4_dim)
    };
};

}
