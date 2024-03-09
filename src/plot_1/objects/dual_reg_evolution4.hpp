///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "plot_1/objects/reg_evolution4.hpp"

namespace Pcr3bpProof
{

class DualRegEvolution
{
public:
    struct Params {};

    DualRegEvolution(
        Lyra::Core3dObjects& core_objects_ref,
        Pcr3bp::SetupParameters<RMap> const & setup,
        Manifold4_Transformation const & transformation_ref,
        RVector initial_point,
        double time,
        size_t point_count,
        float thickness)
    {
        {
            const RegEvolution4::Param param = {
                setup,
                initial_point,
                time,
                point_count,
                transformation_ref,
                thickness,
                true
            };

            m_reg_evolution_pos = std::make_unique<RegEvolution4>(std::ref(core_objects_ref), std::cref(param));
        }
        {
            const RegEvolution4::Param param = {
                setup,
                initial_point,
                time,
                point_count,
                transformation_ref,
                thickness,
                false
            };

            m_reg_evolution_neg = std::make_unique<RegEvolution4>(std::ref(core_objects_ref), std::cref(param));
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
    
private:
    std::unique_ptr<RegEvolution4> m_reg_evolution_pos {};
    std::unique_ptr<RegEvolution4> m_reg_evolution_neg {};
};

}
