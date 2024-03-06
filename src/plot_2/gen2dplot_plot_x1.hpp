///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>

#include "tools/test_tools.hpp"
#include "pcr3bp_obsolete/pcr3bp_reg_poincare_coord.hpp"

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
        m_core_ref.register_surface2(&m_func);
    }

private:
    class Func final : public Lyra::Surface2
    {
    public:
        Func() : Surface2({
            Leo::RulerSet<2>({
                Leo::Ruler<>(-2.5, 2.5, 200, 1),
                Leo::Ruler<>(-2.5, 2.5, 200, 1)
            })
            })
            , m_setup()
            , m_xieta(m_setup)
        {
            Aries::HrClock clock;
            clock.start();

            auto func = [this](const std::array<double, 2>& in) -> std::array<float, 2>
            {
                try
                {
                    const RVector arg = { in[0], in[1] };

                    if ( std::abs(arg.euclNorm() - 2.0) > 1e-6 )
                    {
                        const RVector img = m_xieta(m_xieta( RVector{ in[0], in[1] } ));
                        return { float(img[0]), float(img[1]) };
                    }
                    else
                    {
                        return { NAN, NAN };
                    }
                }
                catch (...)
                {
                    return { NAN, NAN };
                }
            };

            fill(func);

            clock.stop();
            std::cout << "Computation time: " << clock.get_s() << '\n';
        }

    private:
        Pcr3bpSetupValues<RMap> m_setup;
        Pcr3bpRegPoincarePositiveU_XiEta_CE<RMap> m_xieta;

    } m_func;

    Lyra::Core2d& m_core_ref;
};

}
