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


namespace Pcr3bpProof
{

class CoreInterior : CoreInteriorBaseRhez_u_24
{
private:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    std::unique_ptr<RegEvolution4> m_reg_evolution {};
    std::unique_ptr<RegEvolution4> m_reg_evolution_2 {};
    std::unique_ptr<RegEvolution4> m_reg_evolution_3 {};
    std::unique_ptr<RegEvolution4> m_reg_evolution_4 {};

    CoveringRelationsSetup m_covering_relations_setup {};

    std::list<HL_Map> m_origins {};

    std::unique_ptr<SectionPlot4_CE> m_short_path_section_CE {};
    std::unique_ptr<SectionPlot4_CE> m_long_path_section_CE {};

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

        for (auto& hp : m_hset_parameters_list)
        {
            hp.serialize(std::cout);
        }
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

        if (show_periodic_orbit)
        {
            const RVector initial_point = CapdUtils::Concat<MapT>::concat_vectors({ m_basic_objects.m_parameters.get_initial_point(), RVector{ h } });

            reload_reg_evolution(
                m_reg_evolution,
                m_reg_evolution_2,
                initial_point,
                0.908943,
                reg_evo_point_count,
                reg_evo_thickness);
        }
        else
        {
            m_reg_evolution.reset();
            m_reg_evolution_2.reset();
        }

        if (show_homoclinic_orbit)
        {
            const RVector initial_point = { 1.265830729, 0.0, 0.0, 0.1201350685, -0.711058691 };

            reload_reg_evolution(
                m_reg_evolution_3,
                m_reg_evolution_4,
                initial_point,
                2.6362,
                reg_evo_point_count,
                reg_evo_thickness);
        }
        else
        {
            m_reg_evolution_3.reset();
            m_reg_evolution_4.reset();
        }
        
        m_origins.clear();

        using Coordsys = CapdUtils::LocalCoordinateSystem<IMap>;
        auto periodic_orbit_coordsys_vector = m_covering_relations_setup.get_periodic_orbit_coordsys();
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

        auto homoclinic_orbit_coordsys_vector = m_covering_relations_setup.get_homoclinic_orbit_coordsys();
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

        if (centerpoint_index != -1)
        {
            if (centerpoint_index < periodic_orbit_coordsys_vector.size())
            {
                this->set_offset( convert<4>( periodic_orbit_coordsys_vector.at(centerpoint_index).get_origin() ) );
            }
            else if (centerpoint_index < periodic_orbit_coordsys_vector.size() + homoclinic_orbit_coordsys_vector.size())
            {
                size_t const idx = centerpoint_index - periodic_orbit_coordsys_vector.size();
                this->set_offset( convert<4>( periodic_orbit_coordsys_vector.at(idx).get_origin() ) );
            }
        }
        else
        {
            this->set_offset(Lyra::Point4d());
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
    }

    void set_rotation_4d(Leo::Matrix4f const & matrix)
    {
        CoreInteriorBaseRhez_u_24::set_rotation_4d(matrix);


        if (m_reg_evolution)
        {
            m_reg_evolution->refresh();
        }

        if (m_reg_evolution_2)
        {
            m_reg_evolution_2->refresh();
        }

        if (m_reg_evolution_3)
        {
            m_reg_evolution_3->refresh();
        }

        if (m_reg_evolution_4)
        {
            m_reg_evolution_4->refresh();
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
    }

private:
    std::list<CapdUtils::HsetParameters> m_hset_parameters_list {};
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
