///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "hset_parameters_to_coordsys_converter.hpp"

namespace Pcr3bpProof
{
    
HsetParametersToCoordsysConverter::HsetParametersToCoordsysConverter(
    const std::list<CapdUtils::HsetParameters>& hset_parameters_list,
    const std::vector<Coordsys>& periodic_orbit_coordsys_vector,
    const std::vector<Coordsys>& homoclinic_orbit_coordsys_vector)
{
    for (CapdUtils::HsetParameters const & hp : hset_parameters_list)
    {
        m_container[&hp] = find_coordsys(hp, periodic_orbit_coordsys_vector, homoclinic_orbit_coordsys_vector);
    }
}

HsetParametersToCoordsysConverter::Coordsys const & HsetParametersToCoordsysConverter::get_coordsys(CapdUtils::HsetParameters const & hp) const
{
    return *m_container.at(&hp);
}

const HsetParametersToCoordsysConverter::Coordsys* HsetParametersToCoordsysConverter::find_coordsys(
    CapdUtils::HsetParameters const & hp,
    const std::vector<Coordsys>& periodic_orbit_coordsys_vector,
    const std::vector<Coordsys>& homoclinic_orbit_coordsys_vector)
{
    const std::array<double, 4>& coordsys_origin = hp.coordsys_origin;

    const Coordsys* ret = nullptr;
    double distance = INFINITY;

    const RVector co =  { coordsys_origin[0], coordsys_origin[1], coordsys_origin[2], coordsys_origin[3] };

    CapdUtils::MaxNorm<RMap> norm {};

    for (const Coordsys& cs : homoclinic_orbit_coordsys_vector)
    {
        const double d = norm( CapdUtils::vector_cast<RVector>(cs.get_origin()) - co);
        if (d < distance)
        {
            ret = &cs;
            distance = d;
        }
    }

    for (const Coordsys& cs : periodic_orbit_coordsys_vector)
    {
        const double d = norm( CapdUtils::vector_cast<RVector>(cs.get_origin()) - co);
        if (d < distance)
        {
            ret = &cs;
            distance = d;
        }
    }

    return ret;
};

}
