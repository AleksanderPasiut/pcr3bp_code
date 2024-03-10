///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "plot_1/objects/reg_evolution4.hpp"
#include "plot_2/objects/object_with_params.hpp"

namespace Pcr3bpProof
{

class DualRegEvolution
{
public:
    struct Params
    {
        Pcr3bp::SetupParameters<RMap> const & setup;
        RVector initial_point;
        double time;
        size_t point_count;
        float thickness;
    };

    DualRegEvolution(
        Lyra::Core3dObjects& core_objects_ref,
        Manifold4_Transformation const & transformation_ref,
        Params const & params)
    {
        {
            const RegEvolution4::Params param = {
                params.setup,
                params.initial_point,
                params.time,
                params.point_count,
                params.thickness,
                true
            };

            m_reg_evolution_pos = std::make_unique<RegEvolution4>(
                std::ref(core_objects_ref),
                std::cref(transformation_ref),
                std::cref(param));
        }
        {
            const RegEvolution4::Params param = {
                params.setup,
                params.initial_point,
                params.time,
                params.point_count,
                params.thickness,
                false
            };

            m_reg_evolution_neg = std::make_unique<RegEvolution4>(
                std::ref(core_objects_ref),
                std::cref(transformation_ref),
                std::cref(param));
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

class DualRegEvolutionNew
{
public:
    struct Params
    {
        Pcr3bp::SetupParameters<RMap> const & setup;
        RVector initial_point;
        double time;
        size_t point_count;
        float thickness;
    };

    DualRegEvolutionNew(
        Lyra::Core3d& core_ref,
        Manifold4_Transformation const & transformation_ref)
            : m_reg_evolution_pos(core_ref, transformation_ref)
            , m_reg_evolution_neg(core_ref, transformation_ref)
    {}

    void rebuild(Params const & params)
    {
        {
            const RegEvolution4::Params params_internal = {
                params.setup,
                params.initial_point,
                params.time,
                params.point_count,
                params.thickness,
                true
            };

            m_reg_evolution_pos.rebuild(params_internal);
        }
        {
            const RegEvolution4::Params params_internal = {
                params.setup,
                params.initial_point,
                params.time,
                params.point_count,
                params.thickness,
                false
            };

            m_reg_evolution_neg.rebuild(params_internal);
        }
    }

    void highlight(bool arg)
    {
        m_reg_evolution_pos.highlight(arg);
        m_reg_evolution_neg.highlight(arg);
    }        

    void hide()
    {
        m_reg_evolution_pos.hide();
        m_reg_evolution_neg.hide();
    }

    void refresh()
    {
        m_reg_evolution_pos.refresh();
        m_reg_evolution_neg.refresh();
    }

    void heartbeat()
    {
        m_reg_evolution_pos.heartbeat();
        m_reg_evolution_neg.heartbeat();
    }
    
private:
    Renderable4d_WithParams<RegEvolution4> m_reg_evolution_pos;
    Renderable4d_WithParams<RegEvolution4> m_reg_evolution_neg;
};

}
