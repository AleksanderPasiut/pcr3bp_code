///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "covering_relations_setup.hpp"

namespace Pcr3bpProof
{

template<typename MapT>
class CoveringRelationsTestBase
{
protected:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = CapdUtils::LocalCoordinateSystem<MapT>;

    CoveringRelationsTestBase(const CoveringRelationsSetup& setup)
        : m_periodic_orbit_coordsys(setup.get_periodic_orbit_coordsys())
        , m_homoclinic_orbit_coordsys(setup.get_homoclinic_orbit_coordsys())
    {}


    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    const std::vector<Coordsys> m_periodic_orbit_coordsys;
    const std::vector<Coordsys> m_homoclinic_orbit_coordsys;

    const ScalarType m_gain_factor { 75e-11 };
};

}
