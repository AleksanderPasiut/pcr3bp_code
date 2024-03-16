///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "origins_from_coordsys_container.hpp"

namespace Pcr3bpProof
{

OriginsFromCoordsysContainer::OriginsFromCoordsysContainer(std::vector<Coordsys> const & coordsys_vector, Lyra::Core3d& core_ref, Manifold4_Transformation const & transformation_ref)
    : m_coordsys_vector(coordsys_vector)
{
    m_container.reserve( coordsys_vector.size() );

    for (Coordsys const & coordsys : coordsys_vector)
    {
        m_container.emplace_back(std::ref(core_ref), std::cref(transformation_ref));
    }
}

void OriginsFromCoordsysContainer::rebuild(Params const & params)
{
    if (params.thickness != m_params.thickness)
    {
        auto jt = m_coordsys_vector.begin();
        for (auto it = m_container.begin(); it != m_container.end(); ++it, ++jt)
        {
            auto & point_renderable = *it;
            Coordsys const & coordsys = *jt;

            HL_Map::Params const point_params = {
                .U = CapdUtils::vector_cast<RVector>( coordsys.get_origin() ),
                .thickness = params.thickness
            };

            point_renderable.rebuild(point_params);
        }
    }

    if (params.highlight_idx != m_params.highlight_idx)
    {
        m_container[m_params.highlight_idx].highlight(false);
        m_container[params.highlight_idx].highlight(true);
    }

    m_params = params;
}

void OriginsFromCoordsysContainer::refresh()
{
    for (auto& point : m_container)
    {
        point.refresh();
    }
}

void OriginsFromCoordsysContainer::hide()
{
    for (auto& point : m_container)
    {
        point.hide();
    }
}

void OriginsFromCoordsysContainer::heartbeat()
{
    for (auto& point : m_container)
    {
        point.heartbeat();
    }
}

}
