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

#include "rhez_u_24_core_interior_base.hpp"

// #include "plot_1/objects/reg_evolution4.hpp"
#include "plot_1/objects/dual_reg_evolution4.hpp"

#include "plot_2/objects/hl_map.hpp"
#include "plot_2/objects/hset_renderable.hpp"
#include "plot_2/objects/collision_manifold.hpp"

#include "plot_2/objects/object_with_params.hpp"

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

    CapdUtils::HsetParameters periodic_orbit_2_arg_hset
    {
        .type = CapdUtils::HsetType::Argument,
        .coordsys_origin = periodic_hset_parameters_list.rbegin()->coordsys_origin,
        .coordinates = periodic_hset_parameters_list.begin()->coordinates
    };

    std::list<CapdUtils::HsetParameters> hset_parameters_list {};
    hset_parameters_list.splice(hset_parameters_list.end(), homoclinic_hset_parameters_list);
    hset_parameters_list.splice(hset_parameters_list.end(), periodic_hset_parameters_list);
    hset_parameters_list.splice(hset_parameters_list.end(), jump_hset_parameters_list);

    hset_parameters_list.emplace_back(periodic_orbit_2_arg_hset);


    // for (auto& hp : hset_parameters_list)
    // {
    //     hp.serialize(std::cout);
    // }
    return hset_parameters_list;
}



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

        const bool show_collision_manifold = this->get_param(idx++);

        const bool show_periodic_orbit = this->get_param(idx++);
        const bool show_homoclinic_orbit = this->get_param(idx++);

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

        if (show_collision_manifold)
        {
            CollisionManifold::Params const params
            {
                .thickness = reg_evo_thickness
            };

            m_collision_manifold.show(params);
        }
        else
        {
            m_collision_manifold.hide();
        }

        m_dual_reg_evolution_list.clear();

        if (show_periodic_orbit)
        {
            const RVector initial_point = CapdUtils::Concat<MapT>::concat_vectors({ m_basic_objects.m_parameters.get_initial_point(), RVector{ h } });

            DualRegEvolutionNew::Params const params = {
                .setup = m_basic_objects.m_setup,
                .initial_point = initial_point, 
                .time = 0.908943,
                .point_count = reg_evo_point_count,
                .thickness = reg_evo_thickness
            };

            m_periodic_orbit.show(params);
        }
        else
        {
            m_periodic_orbit.hide();
        }

        if (show_homoclinic_orbit)
        {
            const RVector initial_point = { 1.265830729, 0.0, 0.0, 0.1201350685, -0.711058691 };

            DualRegEvolutionNew::Params const params = {
                .setup = m_basic_objects.m_setup,
                .initial_point = initial_point, 
                .time = 2.6362,
                .point_count = reg_evo_point_count,
                .thickness = reg_evo_thickness
            };

            m_homoclinic_orbit.show(params);
        }
        else
        {
            m_homoclinic_orbit.hide();
        }

        if (show_periodic_orbit_local)
        {
            for (Coordsys const& coordsys : periodic_orbit_coordsys_vector)
            {
                DualRegEvolution::Params const params = {
                    .setup = m_basic_objects.m_setup,
                    .initial_point = CapdUtils::vector_cast<RVector>( coordsys.get_origin() ), 
                    .time = evolution_time,
                    .point_count = reg_evo_point_count,
                    .thickness = reg_evo_thickness
                };

                m_dual_reg_evolution_list.emplace_back(
                    std::ref(get_core_ref().get_objects()),
                    std::cref(this->get_transformation()),
                    std::cref(params));
            }
        }

        if (show_homoclinic_orbit_local)
        {
            for (Coordsys const& coordsys : homoclinic_orbit_coordsys_vector)
            {
                DualRegEvolution::Params const params = {
                    .setup = m_basic_objects.m_setup,
                    .initial_point = CapdUtils::vector_cast<RVector>( coordsys.get_origin() ), 
                    .time = evolution_time,
                    .point_count = reg_evo_point_count,
                    .thickness = reg_evo_thickness
                };

                m_dual_reg_evolution_list.emplace_back(
                    std::ref(get_core_ref().get_objects()),
                    std::cref(this->get_transformation()),
                    std::cref(params));
            }
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

        for (auto& evo : m_dual_reg_evolution_list)
        {
            evo.refresh();
        }

        for (auto& ptdbg : m_origins)
        {
            ptdbg.refresh();
        }

        for (auto& hs : m_h_sets)
        {
            hs.refresh();
        }
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

    std::list<DualRegEvolution> m_dual_reg_evolution_list {};
    
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
};

}

int main(int argc, char* argv[])
{
    using Core = Taurus::DefaultCore4d<Pcr3bpProof::CoreInterior>;
    Serpent::SglHostWindow<Core> window(
        Pcr3bpProof::create_window_properties("RHEZ U 24"), argc, argv);

    window.show();
    window.run();
    return 0;
}
