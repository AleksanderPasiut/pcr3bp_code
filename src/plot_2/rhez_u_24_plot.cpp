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

#include "plot_1/objects/reg_evolution4.hpp"
#include "plot_2/objects/hl_map.hpp"
#include "plot_2/objects/section_plot4_ce.hpp"
#include "plot_2/objects/hset_renderable.hpp"


namespace Pcr3bpProof
{

struct DualRegEvolution
{
    std::unique_ptr<RegEvolution4> m_reg_evolution_pos {};
    std::unique_ptr<RegEvolution4> m_reg_evolution_neg {};

    DualRegEvolution(
        Lyra::Core3d& core_ref,
        Pcr3bp::SetupParameters<RMap> setup,
        Manifold4_Transformation const & transformation_ref,
        RVector ret,
        double time,
        size_t point_count,
        float thickness)
    {
        {
            const RegEvolution4::Param param = {
                setup,
                ret,
                time,
                point_count,
                transformation_ref,
                thickness,
                true
            };

            m_reg_evolution_pos = std::make_unique<RegEvolution4>(std::ref(core_ref), std::cref(param));
        }
        {
            const RegEvolution4::Param param = {
                setup,
                ret,
                time,
                point_count,
                transformation_ref,
                thickness,
                false
            };

            m_reg_evolution_neg = std::make_unique<RegEvolution4>(std::ref(core_ref), std::cref(param));
        }
    }

    void refresh()
    {
        if (m_reg_evolution_pos)
        {
            m_reg_evolution_pos->refresh();
        }

        if (m_reg_evolution_neg)
        {
            m_reg_evolution_neg->refresh();
        }
    }
};

class CoreInterior : CoreInteriorBaseRhez_u_24
{
private:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = CapdUtils::LocalCoordinateSystem<IMap>;

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    std::list<DualRegEvolution> m_dual_reg_evolution_list {};
    
    std::unique_ptr<RegEvolution4> m_reg_evolution_3 {};
    std::unique_ptr<RegEvolution4> m_reg_evolution_4 {};

    CoveringRelationsSetup m_covering_relations_setup {};

    std::list<HL_Map> m_origins {};

    std::unique_ptr<SectionPlot4_CE> m_short_path_section_CE {};
    std::unique_ptr<SectionPlot4_CE> m_long_path_section_CE {};

    std::list<HsetRenderable> m_h_sets {};

public:
    CoreInterior(Lyra::Core3d& core_ref) : CoreInteriorBaseRhez_u_24(core_ref)
    {
        auto load_hset_parameters_list = [](std::string file_name) -> std::list<CapdUtils::HsetParameters>
        {
            std::ifstream ifs(file_name);
            if (ifs)
            {
                return CapdUtils::deserialize_hset_parameters_list(ifs);
            }

            throw std::logic_error("Failed to find hset parameters file!");
        };
        auto homoclinic_hset_parameters_list = load_hset_parameters_list("homoclinic_coverings_hset_parameters.csv");
        auto periodic_hset_parameters_list = load_hset_parameters_list("periodic_coverings_hset_parameters.csv");
        auto jump_hset_parameters_list = load_hset_parameters_list("jump_coverings_hset_parameters.csv");

        m_hset_parameters_list = std::list<CapdUtils::HsetParameters>();
        m_hset_parameters_list.splice(m_hset_parameters_list.end(), homoclinic_hset_parameters_list);
        m_hset_parameters_list.splice(m_hset_parameters_list.end(), periodic_hset_parameters_list);
        m_hset_parameters_list.splice(m_hset_parameters_list.end(), jump_hset_parameters_list);

        // for (auto& hp : m_hset_parameters_list)
        // {
        //     hp.serialize(std::cout);
        // }
    }

    void reload_reg_evolution(
        std::unique_ptr<RegEvolution4>& reg_evolution_pos, 
        std::unique_ptr<RegEvolution4>& reg_evolution_neg, 
        RVector ret, double time, size_t point_count, float thickness)
    {
        {
            const RegEvolution4::Param param = {
                m_basic_objects.m_setup,
                ret,
                time,
                point_count,
                this->get_transformation(),
                thickness,
                true
            };

            reg_evolution_pos = std::make_unique<RegEvolution4>(std::ref(get_core_ref()), std::cref(param));
        }
        {
            const RegEvolution4::Param param = {
                m_basic_objects.m_setup,
                ret,
                time,
                point_count,
                this->get_transformation(),
                thickness,
                false
            };

            reg_evolution_neg = std::make_unique<RegEvolution4>(std::ref(get_core_ref()), std::cref(param));
        }
    }

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        CoreInteriorBase::set_param(packet_vector);

        int idx = 0;
        const double point_thickness = this->get_param(idx++);

        const size_t reg_evo_point_count = this->get_param(idx++);
        const float reg_evo_thickness = this->get_param(idx++);

        const int select_short_path_section_CE = this->get_param(idx++);
        const int select_long_path_section_CE = this->get_param(idx++);

        const double section_span = this->get_param(idx++);
        // const double section_scale = this->get_param(idx++);

        const double evolution_time = this->get_param(idx++);

        const double scale = this->get_param(idx++);

        const bool show_periodic_orbit = this->get_param(idx++);
        const bool show_homoclinic_orbit = this->get_param(idx++);

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

        m_dual_reg_evolution_list.clear();

        if (show_periodic_orbit)
        {
            const RVector initial_point = CapdUtils::Concat<MapT>::concat_vectors({ m_basic_objects.m_parameters.get_initial_point(), RVector{ h } });

            m_dual_reg_evolution_list.emplace_back(
                std::ref(get_core_ref()),
                m_basic_objects.m_setup,
                std::cref(this->get_transformation()),
                initial_point,
                0.908943,
                reg_evo_point_count,
                reg_evo_thickness);
        }

        if (show_homoclinic_orbit)
        {
            const RVector initial_point = { 1.265830729, 0.0, 0.0, 0.1201350685, -0.711058691 };

            m_dual_reg_evolution_list.emplace_back(
                std::ref(get_core_ref()),
                m_basic_objects.m_setup,
                std::cref(this->get_transformation()),
                initial_point,
                2.6362,
                reg_evo_point_count,
                reg_evo_thickness);
        }
        
        m_origins.clear();

        
        if (show_periodic_orbit_origins)
        {
            size_t idx = 0;
            for (Coordsys const& coordsys : periodic_orbit_coordsys_vector)
            {
                m_origins.emplace_back(
                    std::ref(get_core_ref()),
                    CapdUtils::vector_cast<RVector>( coordsys.get_origin() ),
                    std::cref( this->get_transformation() ),
                    point_thickness * (idx == highlight_periodic_orbit_point ? 2.0 : 1.0)
                );

                ++idx;
            }
        }

        if (show_homoclinic_orbit_origins)
        {
            size_t idx = 0;
            for (Coordsys const& coordsys : homoclinic_orbit_coordsys_vector)
            {
                m_origins.emplace_back(
                    std::ref(get_core_ref()),
                    CapdUtils::vector_cast<RVector>( coordsys.get_origin() ),
                    std::cref( this->get_transformation() ),
                    point_thickness * (idx == highlight_homoclinic_orbit_point ? 2.0 : 1.0)
                );
                ++idx;
            }
        }

        if (select_short_path_section_CE >= 0)
        {
            const SectionPlot4_CE::Param param
            {
                std::ref(m_basic_objects),
                CapdUtils::LocalCoordinateSystem<MapT>::convert_from( periodic_orbit_coordsys_vector.at(select_short_path_section_CE) ),
                section_span,
                11,
                std::cref(this->get_transformation()),
                reg_evo_thickness
            };

            m_short_path_section_CE = std::make_unique<SectionPlot4_CE>(
                std::ref(get_core_ref()), std::cref(param));
        }
        else
        {
            m_short_path_section_CE.reset();
        }

        if (select_long_path_section_CE >= 0)
        {
            const SectionPlot4_CE::Param param
            {
                std::ref(m_basic_objects),
                CapdUtils::LocalCoordinateSystem<MapT>::convert_from( homoclinic_orbit_coordsys_vector.at(select_long_path_section_CE) ),
                section_span,
                11,
                std::cref(this->get_transformation()),
                reg_evo_thickness
            };

            m_long_path_section_CE = std::make_unique<SectionPlot4_CE>(
                std::ref(get_core_ref()), std::cref(param));
        }
        else
        {
            m_long_path_section_CE.reset();
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

        for (auto& evo : m_dual_reg_evolution_list)
        {
            evo.refresh();
        }

        for (auto& ptdbg : m_origins)
        {
            ptdbg.refresh();
        }

        if (m_short_path_section_CE)
        {
            m_short_path_section_CE->refresh();
        }

        if (m_long_path_section_CE)
        {
            m_long_path_section_CE->refresh();
        }

        for (auto& hs : m_h_sets)
        {
            hs.refresh();
        }
    }

private:
    std::list<CapdUtils::HsetParameters> m_hset_parameters_list {};

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
