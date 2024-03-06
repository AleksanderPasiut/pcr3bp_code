///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>

#include "tools/test_tools.hpp"
#include "pcr3bp_obsolete/pcr3bp_reg_poincare_positive_u_xieta.hpp"

#include <aries/hrclock.hpp>

namespace Ursa
{

class X
{
public:
    using MapT = RMap;
    using ScalarType = MapT::ScalarType;
    using VectorType = MapT::VectorType;
    using MatrixType = MapT::MatrixType;

    X(Lyra::Core2d& core_ref, const std::vector<double>& param)
        : m_core_ref(core_ref)
    {
        m_core_ref.register_surface(&m_func);
        m_core_ref.show_colorbar(true);
    }

private:
    class Func final : public Lyra::Surface
    {
    public:
        Func() : Surface({
            Leo::RulerSet<2>({
                Leo::Ruler<>(-1.5, 1.5, 300, 1),
                Leo::Ruler<>(-1.5, 1.5, 300, 1)
            })
            })
            , m_setup()
            , m_xieta(m_setup)
        {
            Aries::HrClock clock;
            clock.start();

            auto func = [this](const std::array<double, 2>& in) -> std::array<float, 1>
            {
                try
                {
                    const RVector arg = { in[0], in[1] };

                    if ( std::abs(arg.euclNorm() - 2.0) > 1e-6 )
                    {
                        const RVector img = m_xieta(m_xieta( RVector{ in[0], in[1] } )) - arg;
                        const double value = img[0]*img[0] + img[1]*img[1];
                        return { float(value) };
                    }
                    else
                    {
                        return { NAN };
                    }
                }
                catch (...)
                {
                    return { NAN };
                }
            };

            fill(func);

            clock.stop();
            std::cout << "Computation time: " << clock.get_s() << '\n';
        }

    private:
        Pcr3bp::SetupParameters<RMap> m_setup;
        Pcr3bpRegPoincarePositiveU_XiEta_CE<RMap> m_xieta;

    } m_func;

    Lyra::Core2d& m_core_ref;
};

}
