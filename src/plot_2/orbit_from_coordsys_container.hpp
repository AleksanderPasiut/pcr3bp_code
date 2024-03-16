///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <proof/covering_relations_setup.hpp>

#include "plot_1/objects/dual_reg_evolution4.hpp"

#include <vector>

namespace Pcr3bpProof
{

class OrbitFromCoordsysContainer
{
public:
    struct Params
    {
        Pcr3bp::SetupParameters<RMap> const & setup;
        double time;
        size_t point_count;
        float thickness;
        double h;
    };

    using Coordsys = CapdUtils::LocalCoordinateSystem<IMap>;

    OrbitFromCoordsysContainer(std::vector<Coordsys> const & coordsys_vector, Lyra::Core3d& core_ref, Manifold4_Transformation const & transformation_ref);

    void rebuild(Params const & params);

    void refresh();

    void hide();

private:
    std::vector<Coordsys> const & m_coordsys_vector;
    std::vector<DualRegEvolution> m_container {};
};

}
