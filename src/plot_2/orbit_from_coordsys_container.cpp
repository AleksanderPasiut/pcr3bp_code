///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// #include <taurus/default_core_4d.hpp>

// #include <tools/types.hpp>
// #include <tools/logging/h_set_parameters.hpp>

// #include <proof/pcr3bp_reg_basic_objects.hpp>
// #include <proof/covering_relations_setup.hpp>

// #include <capd_utils/concat.hpp>

// #include "plot_common/window_properties.hpp"

// #include "plot_1/objects/dual_reg_evolution4.hpp"


#include "orbit_from_coordsys_container.hpp"

namespace Pcr3bpProof
{

OrbitFromCoordsysContainer::OrbitFromCoordsysContainer(std::vector<Coordsys> const & coordsys_vector, Lyra::Core3d& core_ref, Manifold4_Transformation const & transformation_ref)
    : m_coordsys_vector(coordsys_vector)
{
    m_container.reserve( coordsys_vector.size() );

    for (Coordsys const & coordsys : coordsys_vector)
    {
        m_container.emplace_back(std::ref(core_ref), std::cref(transformation_ref));
    }
}

void OrbitFromCoordsysContainer::rebuild(Params const & params)
{
    auto h = params.h;

    auto get_initial_point_from_coordsys = [h](const Coordsys& coordsys) -> RVector
    {
        return CapdUtils::Concat<RMap>::concat_vectors({
            CapdUtils::vector_cast<RVector>( coordsys.get_origin() ),
            RVector{ h }
        });
    };

    auto jt = m_coordsys_vector.begin();
    for (auto it = m_container.begin(); it != m_container.end(); ++it, ++jt)
    {
        DualRegEvolutionNew & evolution = *it;
        Coordsys const & coordsys = *jt;

        DualRegEvolutionNew::Params const evolution_params = {
            .setup = params.setup,
            .initial_point = get_initial_point_from_coordsys(coordsys),
            .time = params.time,
            .point_count = params.point_count,
            .thickness = params.thickness
        };

        evolution.rebuild(evolution_params);
    }
}

void OrbitFromCoordsysContainer::refresh()
{
    for (DualRegEvolutionNew& evolution : m_container)
    {
        evolution.refresh();
    }
}

void OrbitFromCoordsysContainer::hide()
{
    for (DualRegEvolutionNew& evolution : m_container)
    {
        evolution.hide();
    }
}

}
