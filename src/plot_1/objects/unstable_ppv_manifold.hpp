///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "interp.hpp"

namespace Ursa
{

template<typename PpvType>
class UnstableManifold
{
public:
    UnstableManifold(
        Lyra::Core3d& core_ref,
        const Pcr3bpSetupValues<RMap>& setup,
        const RVector& directory,
        const Leo::Color& color_pos,
        const Leo::Color& color_neg,
        bool full = false)
            : m_interp()
            , m_core_ref(core_ref)
            , m_setup(setup)
    {
        InterpParam param;
        param.setup = setup;
        param.direction = directory;
        param.sub_points = 50;
        param.color = color_pos;
        param.line_thickness = 1.0e-3f;
        param.point_thickness = 0.0e-3f;
        param.show_interpolation_nodes = false;

        const size_t n = 25;

        add_interp(param, 3, 8*n,          0.0, exp(-7.494));
        add_interp(param, 5, 3*n, exp(-14.300), exp(-7.568));

        if (full)
        {
            add_interp(param, 6,  8*n, exp(-13.162), exp(-9.368));
            add_interp(param, 7, 20*n, exp(-13.140), exp(-9.390));
            add_interp(param, 8, 12*n, exp(-13.130), exp(-9.399));
            add_interp(param, 9, 10*n, exp(-13.125), exp(-9.429));
        }

        param.direction = -directory;
        param.color = color_neg;

        add_interp(param, 2, 2*n,          0.0, exp( -7.176));
        add_interp(param, 3, 2*n, exp(-12.200), exp( -7.223));
        add_interp(param, 3, 2*n,          0.0, exp(-12.600));

        if (full)
        {
            add_interp(param, 4, 8*n, exp(-12.100), exp(-10.260));
            add_interp(param, 5, 8*n, exp( -8.001), exp( -7.819));
            add_interp(param, 5, 8*n, exp( -7.816), exp( -7.239));
        }
    }

private:
    void
    add_interp(InterpParam& param, size_t level, size_t points, double min_multiplier, double max_multiplier)
    {
        std::cout << "setting " << m_interp.size() << '\n';

        param.level = level;
        param.points = points;
        param.min_multiplier = min_multiplier;
        param.max_multiplier = max_multiplier;
        m_interp.emplace_back(std::ref(m_core_ref), std::cref(param));
    }

    Lyra::Core3d& m_core_ref;
    const Pcr3bp::SetupParameters<RMap>& m_setup;

public:
    std::list<RenderableInterp<PpvType>> m_interp;
};

}
