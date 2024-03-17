///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "load_hset_parameters_list.hpp"

#include <fstream>

namespace Pcr3bpProof
{

std::list<CapdUtils::HsetParameters> load_hset_parameters_list()
{
    auto load_hset_parameters_list_file = [](std::string file_name) -> std::list<CapdUtils::HsetParameters>
    {
        std::ifstream ifs(file_name);
        if (ifs)
        {
            return CapdUtils::deserialize_hset_parameters_list(ifs);
        }

        throw std::logic_error("Failed to find hset parameters file!");
    };
    auto homoclinic_hset_parameters_list = load_hset_parameters_list_file("homoclinic_coverings_hset_parameters.csv");
    auto periodic_hset_parameters_list = load_hset_parameters_list_file("periodic_coverings_hset_parameters.csv");
    auto jump_hset_parameters_list = load_hset_parameters_list_file("jump_coverings_hset_parameters.csv");
    auto parallelogram_hset_parameters_list = load_hset_parameters_list_file("parallelogram_coverings_hset_parameters.csv");

    CapdUtils::HsetParameters periodic_orbit_2_arg_hset
    {
        .type = CapdUtils::HsetType::Argument,
        .coordsys_origin = periodic_hset_parameters_list.rbegin()->coordsys_origin,
        .coordinates = periodic_hset_parameters_list.begin()->coordinates
    };

    CapdUtils::HsetParameters homoclinic_orbit_K_arg_hset
    {
        .type = CapdUtils::HsetType::Argument,
        .coordsys_origin = homoclinic_hset_parameters_list.rbegin()->coordsys_origin,
        .coordinates = periodic_hset_parameters_list.begin()->coordinates
    };

    std::list<CapdUtils::HsetParameters> hset_parameters_list {};
    hset_parameters_list.splice(hset_parameters_list.end(), homoclinic_hset_parameters_list);
    hset_parameters_list.splice(hset_parameters_list.end(), periodic_hset_parameters_list);
    hset_parameters_list.splice(hset_parameters_list.end(), jump_hset_parameters_list);
    hset_parameters_list.splice(hset_parameters_list.end(), parallelogram_hset_parameters_list);

    hset_parameters_list.emplace_back(periodic_orbit_2_arg_hset);
    hset_parameters_list.emplace_back(homoclinic_orbit_K_arg_hset);



    // for (auto& hp : hset_parameters_list)
    // {
    //     hp.serialize(std::cout);
    // }
    return hset_parameters_list;
}

}
