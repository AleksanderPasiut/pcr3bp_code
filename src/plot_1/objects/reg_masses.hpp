///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>
#include "tools/test_tools.hpp"
#include "capd_renderable.hpp"

#include <pcr3bp_basic/setup_parameters.hpp>
#include <pcr3bp_basic/levi_civita_inverse_coordinate_change.hpp>

namespace Pcr3bpProof
{

class RegMasses
{
private:
    using VectorRenderable = CapdVectorRenderable<RVector, 2>;

    VectorRenderable m_mass_2;
    std::unique_ptr<VectorRenderable> m_mass_1A_ptr;
    std::unique_ptr<VectorRenderable> m_mass_1B_ptr;

public:
    RegMasses(
        Lyra::Core2d& core_ref,
        Pcr3bp::SetupParameters<RMap> setup,
        float point_size = 0.01f,
        Leo::Color color_1 = Leo::Color(0.4, 0.5, 0.8),
        Leo::Color color_2 = Leo::Color(0.8, 0.2, 0.2))
            : m_mass_2(core_ref.get_objects(), { 0.0, 0.0 }, color_2)
            , m_mass_1A_ptr()
            , m_mass_1B_ptr()
    {
        LeviCivitaInverseCoordinateChange<RMap> invc(setup.get_x(2));

        const RVector pos_1A = invc({ setup.get_x(1), 0.0, 0.0, 0.0, 0.0 });
        const RVector pos_1B = -pos_1A;

        m_mass_1A_ptr = std::make_unique<VectorRenderable>(std::ref(core_ref.get_objects()), pos_1A, color_1);
        m_mass_1B_ptr = std::make_unique<VectorRenderable>(std::ref(core_ref.get_objects()), pos_1B, color_1);

        m_mass_2.fill(point_size);
        m_mass_1A_ptr->fill(point_size);
        m_mass_1B_ptr->fill(point_size);
    }
};

}