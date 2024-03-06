///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core2d.hpp>

#include "tools/test_tools.hpp"

#include <complex>

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
                Leo::Ruler<>(-2.0, 3.0, 200, 1),
                Leo::Ruler<>(-2.0, 3.0, 200, 1)
            })
            })
        {
            auto func = [](const std::array<double, 2>& in) -> std::array<double, 2>
            {
                using Complex = std::complex<double>;

                const Complex z( in[0], in[1] );
                // const Complex fz = std::pow(z,3) - Complex(1, 0);

                const Complex I = Complex(0, 1);

                // const Complex fz = (z*z - Complex(1, 0))*pow(z - Complex(2,0) - I, 2) / (z*z + 2.0*Complex(1,1));
                const Complex fz = z;

                return { fz.real(), fz.imag() };
            };

            fill(func);
        }

    } m_func;

    Lyra::Core2d& m_core_ref;
};

}
