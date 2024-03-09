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
    
private:
    std::unique_ptr<RegEvolution4> m_reg_evolution_pos {};
    std::unique_ptr<RegEvolution4> m_reg_evolution_neg {};
};

}
