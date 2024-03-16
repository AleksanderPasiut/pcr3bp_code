///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "rhez_u_24_plot.hpp"

namespace Pcr3bpProof
{

CoreInterior::CoreInterior(Lyra::Core3d& core_ref) : CoreInteriorBaseRhez_u_24(core_ref)
{
    // for (CapdUtils::HsetParameters const & hp : m_hset_parameters_list)
    // {
    //     m_h_sets.emplace_back(
    //         std::ref(get_core_ref()),
    //         std::cref(this->get_transformation()));
    // }
}

void CoreInterior::set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
{
    CoreInteriorBase::set_param(packet_vector);

    int idx = 0;
    const float point_thickness = this->get_param(idx++);

    const size_t reg_evo_point_count = this->get_param(idx++);
    const float reg_evo_thickness = this->get_param(idx++);

    const double evolution_time = this->get_param(idx++);

    const double scale = this->get_param(idx++);

    const unsigned show_collision_manifold = this->get_param(idx++);

    const unsigned show_periodic_orbit = this->get_param(idx++);
    const unsigned show_homoclinic_orbit = this->get_param(idx++);

    const bool show_periodic_orbit_local = this->get_param(idx++);
    const bool show_homoclinic_orbit_local = this->get_param(idx++);

    const bool show_periodic_orbit_origins = this->get_param(idx++);
    const bool show_homoclinic_orbit_origins = this->get_param(idx++);

    const size_t highlight_periodic_orbit_point = this->get_param(idx++);
    const size_t highlight_homoclinic_orbit_point = this->get_param(idx++);

    const bool show_arg_h_sets = this->get_param(idx++);
    const bool show_img_h_sets = this->get_param(idx++);
    const bool show_limg_h_sets = this->get_param(idx++);
    const bool show_rimg_h_sets = this->get_param(idx++);

    const size_t centerpoint_index = this->get_param(idx++);

    this->set_scale(pow(10.0, scale));

    const double h = m_basic_objects.m_parameters.get_energy();

    if (centerpoint_index != -1)
    {
        if (centerpoint_index < m_periodic_orbit_coordsys_vector.size())
        {
            this->set_offset( convert_double<4>( m_periodic_orbit_coordsys_vector.at(centerpoint_index).get_origin() ) );
        }
        else if (centerpoint_index < m_periodic_orbit_coordsys_vector.size() + m_homoclinic_orbit_coordsys_vector.size())
        {
            size_t const idx = centerpoint_index - m_periodic_orbit_coordsys_vector.size();
            this->set_offset( convert_double<4>( m_homoclinic_orbit_coordsys_vector.at(idx).get_origin() ) );
        }
    }
    else
    {
        this->set_offset({});
    }

    if (show_collision_manifold > 0)
    {
        CollisionManifold::Params const params
        {
            .thickness = reg_evo_thickness
        };

        m_collision_manifold.rebuild(params);
        m_collision_manifold.highlight(show_collision_manifold > 1);
    }
    else
    {
        m_collision_manifold.hide();
    }

    if (show_periodic_orbit > 0)
    {
        const RVector initial_point = CapdUtils::Concat<MapT>::concat_vectors({ m_basic_objects.m_parameters.get_initial_point(), RVector{ h } });

        DualRegEvolution::Params const params = {
            .setup = m_basic_objects.m_setup,
            .initial_point = initial_point, 
            .time = 0.908943,
            .point_count = reg_evo_point_count,
            .thickness = reg_evo_thickness
        };

        m_periodic_orbit.rebuild(params);
        m_periodic_orbit.highlight(show_periodic_orbit > 1);
    }
    else
    {
        m_periodic_orbit.hide();
    }

    if (show_homoclinic_orbit > 0)
    {
        const RVector initial_point = { 1.265830729, 0.0, 0.0, 0.1201350685, h };

        DualRegEvolution::Params const params = {
            .setup = m_basic_objects.m_setup,
            .initial_point = initial_point, 
            .time = 2.6362,
            .point_count = reg_evo_point_count,
            .thickness = reg_evo_thickness
        };

        m_homoclinic_orbit.rebuild(params);
        m_homoclinic_orbit.highlight(show_homoclinic_orbit > 1);
    }
    else
    {
        m_homoclinic_orbit.hide();
    }

    if (show_periodic_orbit_local)
    {
        OrbitFromCoordsysContainer::Params params
        {
            .setup = m_basic_objects.m_setup,
            .time = evolution_time,
            .point_count = reg_evo_point_count,
            .thickness = reg_evo_thickness,
            .h = h
        };

        m_periodic_orbit_local.rebuild(params);
    }
    else
    {
        m_periodic_orbit_local.hide();
    }

    if (show_homoclinic_orbit_local)
    {
        OrbitFromCoordsysContainer::Params params
        {
            .setup = m_basic_objects.m_setup,
            .time = evolution_time,
            .point_count = reg_evo_point_count,
            .thickness = reg_evo_thickness,
            .h = h
        };

        m_homoclinic_orbit_local.rebuild(params);
    }
    else
    {
        m_homoclinic_orbit_local.hide();
    }

    if (show_periodic_orbit_origins)
    {
        const OriginsFromCoordsysContainer::Params params
        {
            .thickness = point_thickness,
            .highlight_idx = highlight_periodic_orbit_point
        };

        m_periodic_orbit_origins.rebuild(params);
    }
    else
    {
        m_periodic_orbit_origins.hide();
    }

    if (show_homoclinic_orbit_origins)
    {
        const OriginsFromCoordsysContainer::Params params
        {
            .thickness = point_thickness,
            .highlight_idx = highlight_homoclinic_orbit_point
        };

        m_homoclinic_orbit_origins.rebuild(params);
    }
    else
    {
        m_homoclinic_orbit_origins.hide();
    }

    HsetRenderablesContainer::Params const params
    {
        .basic_objects = m_basic_objects,
        .show_arg_h_sets = show_arg_h_sets,
        .show_img_h_sets = show_img_h_sets,
        .show_limg_h_sets = show_limg_h_sets,
        .show_rimg_h_sets = show_rimg_h_sets,
        .reg_evo_thickness = reg_evo_thickness
    };

    m_hset_renderables_container.update(params);

    // auto it = m_h_sets.begin();
    // for (CapdUtils::HsetParameters const & hp : m_hset_parameters_list)
    // {
    //     bool is_visible = false;
    //     is_visible |= (show_arg_h_sets && hp.type == CapdUtils::HsetType::Argument);
    //     is_visible |= (show_img_h_sets && hp.type == CapdUtils::HsetType::Image);
    //     is_visible |= (show_limg_h_sets && hp.type == CapdUtils::HsetType::LeftImage);
    //     is_visible |= (show_rimg_h_sets && hp.type == CapdUtils::HsetType::RightImage);

    //     auto & h_set = *(it++);
    //     if (is_visible)
    //     {
    //         const Coordsys& coordsys_ref = m_hset_parameter_to_coordsys_converter.get_coordsys(hp);

    //         const HsetRenderable::Params params
    //         {
    //             std::ref(m_basic_objects),
    //             CapdUtils::LocalCoordinateSystem<MapT>::convert_from( coordsys_ref ),
    //             hp.coordinates,
    //             3,
    //             5,
    //             reg_evo_thickness
    //         };

    //         h_set.rebuild(std::cref(params));
    //     }
    //     else
    //     {
    //         h_set.hide();
    //     }
    // }
}

void CoreInterior::set_rotation_4d(Leo::Matrix4f const & matrix)
{
    CoreInteriorBaseRhez_u_24::set_rotation_4d(matrix);

    m_collision_manifold.refresh();

    m_periodic_orbit.refresh();
    m_homoclinic_orbit.refresh();

    m_periodic_orbit_origins.refresh();
    m_homoclinic_orbit_origins.refresh();

    m_periodic_orbit_local.refresh();
    m_homoclinic_orbit_local.refresh();

    m_hset_renderables_container.refresh();

    // for (auto& hs : m_h_sets)
    // {
    //     hs.refresh();
    // }
}

void CoreInterior::heartbeat()
{
    m_collision_manifold.heartbeat();

    m_periodic_orbit.heartbeat();
    m_homoclinic_orbit.heartbeat();

    m_periodic_orbit_origins.heartbeat();
    m_homoclinic_orbit_origins.heartbeat();
}

}
