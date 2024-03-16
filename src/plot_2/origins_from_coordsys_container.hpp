///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <proof/covering_relations_setup.hpp>

#include "plot_2/objects/hl_map.hpp"
#include "plot_2/objects/renderable4d_with_params.hpp"

#include <vector>

namespace Pcr3bpProof
{

class OriginsFromCoordsysContainer
{
public:
    struct Params
    {
        float thickness;
        size_t highlight_idx;
    };

    using Coordsys = CapdUtils::LocalCoordinateSystem<IMap>;

    OriginsFromCoordsysContainer(std::vector<Coordsys> const & coordsys_vector, Lyra::Core3d& core_ref, Manifold4_Transformation const & transformation_ref);

    void rebuild(Params const & params);

    void refresh();

    void hide();

private:
    std::vector<Coordsys> const & m_coordsys_vector;
    std::vector<Renderable4d_WithParams<HL_Map>> m_container {};
};

}
