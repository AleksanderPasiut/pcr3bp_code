///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "rhez_u_24_core_interior_base.hpp"

#include "plot_1/objects/dual_reg_evolution4.hpp"
#include "plot_2/objects/hset_renderable.hpp"
#include "plot_2/objects/collision_manifold.hpp"
#include "plot_2/objects/renderable4d_with_params.hpp"

#include "load_hset_parameters_list.hpp"
#include "orbit_from_coordsys_container.hpp"
#include "origins_from_coordsys_container.hpp"
#include "hset_parameters_to_coordsys_converter.hpp"

namespace Pcr3bpProof
{

class CoreInterior : CoreInteriorBaseRhez_u_24
{
public:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = CapdUtils::LocalCoordinateSystem<IMap>;

    CoreInterior(Lyra::Core3d& core_ref);

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector);

    void set_rotation_4d(Leo::Matrix4f const & matrix);

    void heartbeat();

private:
    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    DualRegEvolution m_periodic_orbit
    {
        get_core_ref(),
        this->get_transformation()
    };

    DualRegEvolution m_homoclinic_orbit
    {
        get_core_ref(),
        this->get_transformation()
    };

    Renderable4d_WithParams<CollisionManifold> m_collision_manifold
    {
        get_core_ref(),
        this->get_transformation()
    };

    CoveringRelationsSetup m_covering_relations_setup {};

    std::list<HsetRenderable> m_h_sets {};

    std::list<CapdUtils::HsetParameters> m_hset_parameters_list
    {
        load_hset_parameters_list()
    };

    std::vector<Coordsys> periodic_orbit_coordsys_vector
    {
        m_covering_relations_setup.get_periodic_orbit_coordsys()
    };

    std::vector<Coordsys> homoclinic_orbit_coordsys_vector
    {
        m_covering_relations_setup.get_homoclinic_orbit_coordsys()
    };

    HsetParametersToCoordsysConverter hset_parameter_to_coordsys_converter
    {
        m_hset_parameters_list,
        periodic_orbit_coordsys_vector,
        homoclinic_orbit_coordsys_vector
    };

    OriginsFromCoordsysContainer m_periodic_orbit_origins
    {
        periodic_orbit_coordsys_vector,
        get_core_ref(),
        this->get_transformation()
    };

    OriginsFromCoordsysContainer m_homoclinic_orbit_origins
    {
        homoclinic_orbit_coordsys_vector,
        get_core_ref(),
        this->get_transformation()
    };

    OrbitFromCoordsysContainer m_periodic_orbit_local
    {
        periodic_orbit_coordsys_vector,
        get_core_ref(),
        this->get_transformation()
    };

    OrbitFromCoordsysContainer m_homoclinic_orbit_local
    {
        homoclinic_orbit_coordsys_vector,
        get_core_ref(),
        this->get_transformation()
    };
};

}
