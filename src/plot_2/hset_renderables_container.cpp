///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "hset_renderables_container.hpp"

namespace Pcr3bpProof
{

HsetRenderablesContainer::HsetRenderablesContainer(
    Lyra::Core3d& core_ref,
    const Manifold4_Transformation & transformation_ref,
    std::vector<Coordsys> const & periodic_orbit_coordsys_vector,
    std::vector<Coordsys> const & homoclinic_orbit_coordsys_vector)
        : m_hset_parameter_to_coordsys_converter(
            m_hset_parameters_list,
            periodic_orbit_coordsys_vector,
            homoclinic_orbit_coordsys_vector)
        , m_parallelogram_covering(core_ref, transformation_ref)
{
    m_container.reserve(m_hset_parameters_list.size());
    for (CapdUtils::HsetParameters const & hp : m_hset_parameters_list)
    {
        m_container.emplace_back(
            std::ref(core_ref),
            std::cref(transformation_ref));
    }

    m_parallelogram_coordsys = periodic_orbit_coordsys_vector.at(3);
}

void HsetRenderablesContainer::update(Params const & params)
{
    auto it = m_container.begin();
    for (CapdUtils::HsetParameters const & hp : m_hset_parameters_list)
    {
        bool is_visible = false;
        is_visible |= (params.show_arg_h_sets && hp.type == CapdUtils::HsetType::Argument);
        is_visible |= (params.show_img_h_sets && hp.type == CapdUtils::HsetType::Image);
        is_visible |= (params.show_limg_h_sets && hp.type == CapdUtils::HsetType::LeftImage);
        is_visible |= (params.show_rimg_h_sets && hp.type == CapdUtils::HsetType::RightImage);

        auto & h_set = *(it++);
        if (is_visible)
        {
            const Coordsys& coordsys_ref = m_hset_parameter_to_coordsys_converter.get_coordsys(hp);

            const HsetRenderable::Params params_internal
            {
                params.basic_objects,
                CapdUtils::LocalCoordinateSystem<RMap>::convert_from( coordsys_ref ),
                hp.coordinates,
                3,
                5,
                params.reg_evo_thickness
            };

            h_set.rebuild(std::cref(params_internal));
        }
        else
        {
            h_set.hide();
        }
    }

    if (params.show_parallelogram_h_sets)
    {
        const HsetRenderableParallelogram::Params params_internal
        {
            params.basic_objects,
            CapdUtils::LocalCoordinateSystem<RMap>::convert_from( m_parallelogram_coordsys ),
            { -1.0, 1.0, -1.0, 1.0 },
            3,
            5,
            params.reg_evo_thickness
        };

        m_parallelogram_covering.rebuild(params_internal);
    }
}

void HsetRenderablesContainer::refresh()
{
    for (auto& hs : m_container)
    {
        hs.refresh();
    }

    m_parallelogram_covering.refresh();
}

}
