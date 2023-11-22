///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/newton_method/newton_method.hpp>
#include <carina/implicit_function_differential_equation.hpp>

#include "pcr3bp_std_lyapunov_orbit_lookup.hpp"

namespace Ursa
{

class Pcr3bpStdLyapunovOrbitImplicit
{
public:
    static IVector calculate_py(Real x0, size_t steps, Real& return_time)
    {
        const Pcr3bp::SetupParameters<IMap> setup;

        const IVector PY0 = calculate_initial_py(x0, steps, setup);

        Pcr3bpStdLyapunovOrbitPy<IMap> orbit(setup, x0);
        Carina::NewtonMethod newton(orbit, PY0, 20);
        const IVector PY = newton.get_root();

        return_time = Carina::scalar_cast<Real>( orbit.get_return_time(PY) );
        return PY;
    }

private:
    static IVector calculate_initial_py(Real x0, size_t steps, const Pcr3bp::SetupParameters<IMap>& setup)
    {
        Pcr3bpStdLyapunovOrbitX0Py<IMap> orbit_2(setup);

        Carina::ImplicitFunctionDifferentialEquation<decltype(orbit_2), 1> imp_fun(orbit_2);

        const Pcr3bpStdLyapunovOrbitLookupTable::Entry base_entry
            = Pcr3bpStdLyapunovOrbitLookupTable::get().get_base(x0);

        const IVector PY0 = { imp_fun.evolve(base_entry.x, { base_entry.p }, x0, steps) [0] };
        return PY0;
    }
};

}
