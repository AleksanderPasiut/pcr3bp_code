///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_4d.hpp>

#include <tools/types.hpp>
#include <tools/logging/h_set_parameters.hpp>

#include <proof/pcr3bp_reg_basic_objects.hpp>
#include <proof/covering_relations_setup.hpp>

#include <capd_utils/concat.hpp>

#include "plot_common/window_properties.hpp"
#include "plot_common/timelevel_divisor.hpp"

#include "rhez_u_24_core_interior_base.hpp"

// #include "plot_1/objects/reg_evolution4.hpp"
#include "plot_1/objects/dual_reg_evolution4.hpp"

#include "plot_2/objects/hl_map.hpp"
#include "plot_2/objects/hset_renderable.hpp"
#include "plot_2/objects/collision_manifold.hpp"

#include "plot_2/objects/object_with_params.hpp"
#include "load_hset_parameters_list.hpp"
#include "orbit_from_coordsys_container.hpp"

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

    CoreInterior(Lyra::Core3d& core_ref) : CoreInteriorBaseRhez_u_24(core_ref)
    {}

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
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
            if (centerpoint_index < periodic_orbit_coordsys_vector.size())
            {
                this->set_offset( convert_double<4>( periodic_orbit_coordsys_vector.at(centerpoint_index).get_origin() ) );
            }
            else if (centerpoint_index < periodic_orbit_coordsys_vector.size() + homoclinic_orbit_coordsys_vector.size())
            {
                size_t const idx = centerpoint_index - periodic_orbit_coordsys_vector.size();
                this->set_offset( convert_double<4>( homoclinic_orbit_coordsys_vector.at(idx).get_origin() ) );
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

            DualRegEvolutionNew::Params const params = {
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

            DualRegEvolutionNew::Params const params = {
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

        auto get_initial_point_from_coordsys = [h](const Coordsys& coordsys) -> RVector
        {
            return CapdUtils::Concat<MapT>::concat_vectors({
                CapdUtils::vector_cast<RVector>( coordsys.get_origin() ),
                RVector{ h }
            });
        };

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
        
        m_origins.clear();



        
        if (show_periodic_orbit_origins)
        {
            size_t idx = 0;
            for (Coordsys const& coordsys : periodic_orbit_coordsys_vector)
            {
                HL_Map::Param const param = 
                {
                    .U = CapdUtils::vector_cast<RVector>( coordsys.get_origin() ),
                    .thickness = point_thickness * (idx == highlight_periodic_orbit_point ? 2.0f : 1.0f)
                };

                m_origins.emplace_back(
                    std::ref(get_core_ref().get_objects()),
                    std::cref(this->get_transformation()),
                    std::cref(param) );

                ++idx;
            }
        }

        if (show_homoclinic_orbit_origins)
        {
            size_t idx = 0;
            for (Coordsys const& coordsys : homoclinic_orbit_coordsys_vector)
            {
                HL_Map::Param const param = 
                {
                    .U = CapdUtils::vector_cast<RVector>( coordsys.get_origin() ),
                    .thickness = point_thickness * (idx == highlight_periodic_orbit_point ? 2.0f : 1.0f)
                };

                m_origins.emplace_back(
                    std::ref(get_core_ref().get_objects()),
                    std::cref(this->get_transformation()),
                    std::cref(param) );

                ++idx;
            }
        }

        m_h_sets.clear();

        for (CapdUtils::HsetParameters const & hp : m_hset_parameters_list)
        {
            bool is_visible = false;
            is_visible |= (show_arg_h_sets && hp.type == CapdUtils::HsetType::Argument);
            is_visible |= (show_img_h_sets && hp.type == CapdUtils::HsetType::Image);
            is_visible |= (show_limg_h_sets && hp.type == CapdUtils::HsetType::LeftImage);
            is_visible |= (show_rimg_h_sets && hp.type == CapdUtils::HsetType::RightImage);

            if (is_visible)
            {
                auto find_coordsys = [&](std::array<double, 4> coordsys_origin) -> CapdUtils::LocalCoordinateSystem<IMap>*
                {
                    CapdUtils::LocalCoordinateSystem<IMap>* ret = nullptr;
                    double distance = INFINITY;

                    RVector co =  { coordsys_origin[0], coordsys_origin[1], coordsys_origin[2], coordsys_origin[3] };

                    CapdUtils::MaxNorm<RMap> norm {};

                    for (auto& cs : homoclinic_orbit_coordsys_vector)
                    {
                        double d = norm( CapdUtils::vector_cast<RVector>(cs.get_origin()) - co);
                        if (d < distance)
                        {
                            ret = &cs;
                            distance = d;
                        }
                    }

                    for (auto& cs : periodic_orbit_coordsys_vector)
                    {
                        double d = norm( CapdUtils::vector_cast<RVector>(cs.get_origin()) - co);
                        if (d < distance)
                        {
                            ret = &cs;
                            distance = d;
                        }
                    }

                    return ret;
                };

                CapdUtils::LocalCoordinateSystem<IMap>* coordsys_ptr = find_coordsys(hp.coordsys_origin);

                const HsetRenderable::Param param
                {
                    std::ref(m_basic_objects),
                    CapdUtils::LocalCoordinateSystem<MapT>::convert_from( *coordsys_ptr ),
                    hp.coordinates,
                    3,
                    5,
                    std::cref(this->get_transformation()),
                    reg_evo_thickness
                };

                m_h_sets.emplace_back(
                    std::ref(get_core_ref()),
                    std::cref(param));
            }
        }
    }

    void set_rotation_4d(Leo::Matrix4f const & matrix)
    {
        CoreInteriorBaseRhez_u_24::set_rotation_4d(matrix);

        m_collision_manifold.refresh();

        m_periodic_orbit.refresh();
        m_homoclinic_orbit.refresh();

        m_periodic_orbit_local.refresh();
        m_homoclinic_orbit_local.refresh();

        for (auto& ptdbg : m_origins)
        {
            ptdbg.refresh();
        }

        for (auto& hs : m_h_sets)
        {
            hs.refresh();
        }
    }

    void heartbeat()
    {
        m_collision_manifold.heartbeat();
        m_periodic_orbit.heartbeat();
        m_homoclinic_orbit.heartbeat();
    }

private:
    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    DualRegEvolutionNew m_periodic_orbit
    {
        get_core_ref(),
        this->get_transformation()
    };

    DualRegEvolutionNew m_homoclinic_orbit
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

    std::list<HL_Map> m_origins {};

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

    TimelevelDivisor m_timelevel_divisor {};
};

}

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore4d<Pcr3bpProof::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Pcr3bpProof::create_window_properties("RHEZ U 24"),
        argc,
        argv,
        Leo::Color::White(),
        Serpent::Grid3d::Properties{
            .ruler_x = Leo::Ruler<>(-2, 2, 9, 1),
            .ruler_y = Leo::Ruler<>(-2, 2, 9, 1),
            .ruler_z = Leo::Ruler<>(-2, 2, 9, 1),
            .grid_color = Leo::Color::DarkGray(),
            .grid_width = 0.005f,
            .label_x = L".",
            .label_y = L".",
            .label_z = L".",
            .label_properties{
                .back_color = Leo::Color::White(),
                .fore_color = Leo::Color::DarkGray(),
                .halfheight = 0.05f
            }
        },
        Serpent::Colortray::get_properties_1d(
            Leo::Color::Profile::RedPurpleBlue11,
            Leo::Color::White(),
            Leo::Color::DarkGray(),
            Lyra::Colortray::get_default_placement_1d(),
            Leo::Ruler<double>(-1.0, +1.0, 3, 1))
    );

    window.show();
    window.run();
    return 0;
}
