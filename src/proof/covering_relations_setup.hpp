///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "periodic_orbit_coordsys_generator.hpp"
#include "homoclinic_orbit_origins_initial.hpp"
#include "homoclinic_orbit_origins_generator.hpp"
#include "homoclinic_orbit_coordsys_generator.hpp"

#include "tools/coordsys_utilities.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Generate and store local coordsys for periodic orbit and homoclinic orbit for later use in covering relations
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class CoveringRelationsSetup
{
public:
    std::vector<Carina::LocalCoordinateSystem<IMap>> get_periodic_orbit_coordsys() const
    {
        #if 1

        RegLyapunovCollisionOrbitParameters<IMap> m_parameters {};

        std::vector<Carina::LocalCoordinateSystem<IMap>> periodic_orbit_coordsys(4);

        periodic_orbit_coordsys.at(0) = Carina::LocalCoordinateSystem<IMap>(
            m_parameters.get_initial_point(),
            Carina::matrix_cast<IMatrix>(m_periodic_orbit_coordsys_approx.at(0).get_directions_matrix()) );
        
        periodic_orbit_coordsys.at(1) = Carina::LocalCoordinateSystem<IMap>(
            m_parameters.get_intermediate_point(),
            Carina::matrix_cast<IMatrix>(m_periodic_orbit_coordsys_approx.at(1).get_directions_matrix()) );
        
        periodic_orbit_coordsys.at(2) = Carina::LocalCoordinateSystem<IMap>(
            m_parameters.get_image_point(),
            Carina::matrix_cast<IMatrix>(m_periodic_orbit_coordsys_approx.at(2).get_directions_matrix()) );
        
        periodic_orbit_coordsys.at(3) = Carina::LocalCoordinateSystem<IMap>(
            m_parameters.get_intermediate_point_neg(),
            Carina::matrix_cast<IMatrix>(m_periodic_orbit_coordsys_approx.at(3).get_directions_matrix()) );

        return periodic_orbit_coordsys;

        #else

            return CoordsysVec<IMap>::convert(m_periodic_orbit_coordsys_approx);

        #endif
    }

    std::vector<Carina::LocalCoordinateSystem<IMap>> get_homoclinic_orbit_coordsys() const
    {
        return CoordsysVec<IMap>::convert( m_homoclinic_orbit_coordsys_generator.get_coordsys_container() );
    }

private:
    PeriodicOrbitCoordsysGenerator<RMap> m_periodic_orbit_coordsys_generator_approx {};

    std::vector<Carina::LocalCoordinateSystem<RMap>> m_periodic_orbit_coordsys_approx
    {
        m_periodic_orbit_coordsys_generator_approx.get_coordsys_container()
    };

    HomoclinicOrbitOriginsInitial<RMap> m_homoclinic_orbit_origins_initial {};

    HomoclinicOrbitOriginsGenerator<RMap> m_homoclinic_orbit_origins_generator
    {
        m_homoclinic_orbit_origins_initial
    };

    HomoclinicOrbitCoordsysGenerator<RMap> m_homoclinic_orbit_coordsys_generator
    {
        m_periodic_orbit_coordsys_approx,
        m_homoclinic_orbit_origins_generator.get_points(),
        m_homoclinic_orbit_origins_generator.get_total_expansion_factor_pos(),
        m_homoclinic_orbit_origins_generator.get_total_expansion_factor_neg()
    };
};

}
