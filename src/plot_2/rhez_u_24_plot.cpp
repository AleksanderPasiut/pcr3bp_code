///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <serpent/sgl_host_window.hpp>
#include <taurus/default_core_4d.hpp>

#include <tools/types.hpp>

#if 0

#include <capd_utils/concat.hpp>

#endif

#include "plot_common/window_properties.hpp"

#include "rhez_u_24_core_interior_base.hpp"

#include "plot_1/objects/reg_evolution4.hpp"

#if 0

#include "plot_2/objects/hl_map.hpp"
#include "plot_2/objects/section_plot4.hpp"
#include "plot_2/objects/section_plot4_ce.hpp"

#include <beta/periodic_orbit_parameters.hpp>
#include <beta/covering_relations_setup.hpp>

#endif

namespace Pcr3bpProof
{

class CoreInterior : CoreInteriorBaseRhez_u_24
{
private:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

#if 0

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    std::unique_ptr<RegEvolution4> m_reg_evolution {};
    std::unique_ptr<RegEvolution4> m_reg_evolution_2 {};
    std::array<std::unique_ptr<HL_Map>, 4> m_ptdbg {};

    CoveringRelationsSetup m_covering_relations_setup {};

    std::unique_ptr<SectionPlot4> m_short_path_section {};
    std::unique_ptr<SectionPlot4> m_long_path_section {};

    std::unique_ptr<SectionPlot4_CE> m_short_path_section_CE {};
    std::unique_ptr<SectionPlot4_CE> m_long_path_section_CE {};

#endif

public:
    CoreInterior(Lyra::Core3d& core_ref) : CoreInteriorBaseRhez_u_24(core_ref)
    {}

#if 0

    void reload_reg_evolution(RVector ret, double time, size_t point_count, float thickness)
    {
        const RegEvolution4::Param param = {
            m_basic_objects.m_setup,
            ret,
            time,
            point_count,
            this->get_rotation(),
            thickness,
            true
        };

        m_reg_evolution = std::make_unique<RegEvolution4>(std::ref(get_core_ref()), std::cref(param));
    }

    void reload_reg_evolution_2(RVector ret, double time, size_t point_count, float thickness)
    {
        const RegEvolution4::Param param = {
            m_basic_objects.m_setup,
            ret,
            time,
            point_count,
            this->get_rotation(),
            thickness,
            false
        };

        m_reg_evolution_2 = std::make_unique<RegEvolution4>(std::ref(get_core_ref()), std::cref(param));
    }

#endif

    void set_param(const std::vector<Aquila::ParamPacket<double>>& packet_vector)
    {
        CoreInteriorBase::set_param(packet_vector);

#if 0
        int idx = 0;
        const size_t reg_evo_select = this->get_param(idx++);
        const size_t reg_evo_point_count = this->get_param(idx++);
        const double reg_evo_thickness = this->get_param(idx++);

        const float thickness = this->get_param(idx++);

        const bool show_ghv = this->get_param(idx++) != 0.0;
        // const bool show_mid = this->get_param(idx++) != 0.0;
        // const bool show_neg = this->get_param(idx++) != 0.0;

        const double lyapunov_u0_param = this->get_param(idx++);
        const size_t lyapunov_adj_steps = this->get_param(idx++);

        const int select_short_path_section = this->get_param(idx++);
        const int select_long_path_section = this->get_param(idx++);

        const int select_short_path_section_CE = this->get_param(idx++);
        const int select_long_path_section_CE = this->get_param(idx++);

        const double section_span = this->get_param(idx++);
        // const double section_scale = this->get_param(idx++);

        const double evolution_time = this->get_param(idx++);

        const double h = m_basic_objects.m_parameters.get_energy();

        if (reg_evo_select == -1)
        {
            m_reg_evolution.reset();
            m_reg_evolution_2.reset();
        }

        if (reg_evo_select == 0)
        {
            const RVector initial_point = CapdUtils::Concat<MapT>::concat_vectors({ m_basic_objects.m_parameters.get_initial_point(), RVector{ h } });

            reload_reg_evolution(
                initial_point,
                evolution_time,
                reg_evo_point_count,
                reg_evo_thickness);

            reload_reg_evolution_2(
                initial_point,
                evolution_time,
                reg_evo_point_count,
                reg_evo_thickness);
        }

        if (reg_evo_select == 1)
        {
            const RVector initial_point = { 1.265830729, 0.0, 0.0, 0.1201350685, -0.711058691 };
            const RVector image_point = m_basic_objects.m_parameters.get_image_point_approx();

            reload_reg_evolution(
                initial_point,
                evolution_time,
                reg_evo_point_count,
                reg_evo_thickness);

            reload_reg_evolution_2(
                initial_point,
                evolution_time,
                reg_evo_point_count,
                reg_evo_thickness);
        }

        reload_pos_manifold(h, show_ghv, thickness);
        reload_mid_manifold(h, show_ghv, thickness);
        reload_neg_manifold(h, show_ghv, thickness);

        m_ptdbg[0] = std::make_unique<HL_Map>(
            std::ref(get_core_ref()), m_basic_objects.m_parameters.get_initial_point(), std::cref(this->get_rotation()));

        m_ptdbg[1] = std::make_unique<HL_Map>(
            std::ref(get_core_ref()), m_basic_objects.m_parameters.get_image_point_approx(), std::cref(this->get_rotation()));

        if (select_short_path_section >= 0)
        {
            const SectionPlot4::Param param
            {
                m_covering_relations_setup.get_periodic_orbit_coordsys().at(select_short_path_section),
                1.0,
                std::cref(this->get_rotation()),
                thickness
            };

            m_short_path_section = std::make_unique<SectionPlot4>(
                std::ref(get_core_ref()), std::cref(param));
        }
        else
        {
            m_short_path_section.reset();
        }

        if (select_long_path_section >= 0)
        {
            const SectionPlot4::Param param
            {
                m_covering_relations_setup.get_homoclinic_orbit_coordsys().at(select_long_path_section),
                section_span,
                std::cref(this->get_rotation()),
                thickness
            };

            m_long_path_section = std::make_unique<SectionPlot4>(
                std::ref(get_core_ref()), std::cref(param));
        }
        else
        {
            m_long_path_section.reset();
        }

        if (select_short_path_section_CE >= 0)
        {
            const SectionPlot4_CE::Param param
            {
                std::ref(m_basic_objects),
                m_covering_relations_setup.get_periodic_orbit_coordsys().at(select_short_path_section_CE),
                section_span,
                1.0 / section_span,
                std::cref(this->get_rotation()),
                thickness
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
                m_covering_relations_setup.get_homoclinic_orbit_coordsys().at(select_long_path_section_CE),
                section_span,
                1.0 / section_span,
                std::cref(this->get_rotation()),
                thickness
            };

            m_long_path_section_CE = std::make_unique<SectionPlot4_CE>(
                std::ref(get_core_ref()), std::cref(param));
        }
        else
        {
            m_long_path_section_CE.reset();
        }

#endif
    }

    void set_rotation_4d(Leo::Matrix4f const & matrix)
    {
        CoreInteriorBaseRhez_u_24::set_rotation_4d(matrix);

#if 0

        if (m_reg_evolution)
        {
            m_reg_evolution->refresh();
        }

        if (m_reg_evolution_2)
        {
            m_reg_evolution_2->refresh();
        }

        for (auto& ptdbg : m_ptdbg)
        {
            if (ptdbg)
            {
                ptdbg->refresh();
            }
        }

        if (m_short_path_section)
        {
            m_short_path_section->refresh();
        }

        if (m_long_path_section)
        {
            m_long_path_section->refresh();
        }

        if (m_short_path_section_CE)
        {
            m_short_path_section_CE->refresh();
        }

        if (m_long_path_section_CE)
        {
            m_long_path_section_CE->refresh();
        }

#endif
    }
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
