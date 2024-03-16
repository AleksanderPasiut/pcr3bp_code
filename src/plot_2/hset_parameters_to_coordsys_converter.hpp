///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <tools/logging/h_set_parameters.hpp>
#include <proof/covering_relations_setup.hpp>

#include <map>

namespace Pcr3bpProof
{

class HsetParametersToCoordsysConverter
{
public:
    using Coordsys = CapdUtils::LocalCoordinateSystem<IMap>;
    
    HsetParametersToCoordsysConverter(
        const std::list<CapdUtils::HsetParameters>& hset_parameters_list,
        const std::vector<Coordsys>& periodic_orbit_coordsys_vector,
        const std::vector<Coordsys>& homoclinic_orbit_coordsys_vector);

    Coordsys const & get_coordsys(CapdUtils::HsetParameters const & hp) const;

private:
    const Coordsys* find_coordsys(
        CapdUtils::HsetParameters const & hp,
        const std::vector<Coordsys>& periodic_orbit_coordsys_vector,
        const std::vector<Coordsys>& homoclinic_orbit_coordsys_vector);

    std::map<CapdUtils::HsetParameters const *, Coordsys const *> m_container;
};

}
