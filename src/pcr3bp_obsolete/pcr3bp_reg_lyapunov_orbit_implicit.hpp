///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/newton_method/newton_method.hpp>
#include <carina/implicit_function_differential_equation.hpp>

#include "pcr3bp_reg_lyapunov_orbit_lookup.hpp"

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! (u0, v0, pu0, pv0, u1, v1, pu1, pv1, u2, v2, pu2, pv2, h0) -> (timemap(U0) - U1, poincare(U1) - U2, Gamma(U2))
//! where Gamma is regularized hamiltonian
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LyapunovOrbitRegBaseExtended : Carina::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LyapunovOrbitRegBaseExtended(const Pcr3bp::SetupParameters<MapT>& setup, Real t)
        : m_base(setup, t)
        , m_gamma( Pcr3bp::RegularizedSystem<MapT>::createHamiltonian(2, setup) )
        , m_base_pne(13,
            Carina::IdxList<size_t>::create(0, 13),
            Carina::IdxList<int>::create(11, 0, 10),
            std::ref(m_base))
        , m_gamma_pne(13,
            Carina::IdxList<size_t>{ 4, 5, 6, 7, 12 },
            Carina::IdxList<int>::create(11, 10, 1),
            std::ref(m_gamma))
    {}

    VectorType operator() (const VectorType& vec, MatrixType& mat)
    {
        this->assert_vector_size(vec, 13, "LyapunovOrbitRegBaseExtended vec vector size mismatch!");

        MatrixType m1(11, 13);
        const VectorType v1 = m_base_pne(vec, m1);
        this->assert_matrix_size(m1, 11, 13, "LyapunovOrbitRegBaseExtended m1 matrix size mismatch!");

        MatrixType m2(11, 13);
        const VectorType v2 = m_gamma_pne(vec, m2);
        this->assert_matrix_size(m2, 11, 13, "LyapunovOrbitRegBaseExtended m2 matrix size mismatch!");

        mat = m1 + m2;
        return v1 + v2;
    }

    VectorType operator() (const VectorType& vec)
    {
        MatrixType dummy;
        return (*this)(vec, dummy);
    }

    unsigned dimension() const noexcept
    {
        return 13;
    }

    unsigned imageDimension() const noexcept
    {
        return 11;
    }

private:
    LyapunovOrbitRegBase<MapT> m_base;
    MapT m_gamma;

    Carina::PNE<MapT, decltype(m_base)&> m_base_pne;
    Carina::PNE<MapT, decltype(m_gamma)&> m_gamma_pne;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! TODO
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LyapunovOrbitReg : public Carina::ENP<MapT, LyapunovOrbitRegBaseExtended<MapT>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LyapunovOrbitReg(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType u0, Real t = 1.0)
        : Carina::ENP<MapT, LyapunovOrbitRegBaseExtended<MapT>>(
            VectorType{ u0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
            { -1, -1, -1, 0,
              1,  2,  3,  4,
              -1, -1, -1, -1, 5 },
            { 0, 1, 2, 3, 7, 10 },

            std::cref(setup), t)
    {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! TODO
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LyapunovOrbitReg_Implicit : public Carina::ENP<MapT, LyapunovOrbitRegBaseExtended<MapT>>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LyapunovOrbitReg_Implicit(const Pcr3bp::SetupParameters<MapT>& setup, Real t = 1.0)
        : Carina::ENP<MapT, LyapunovOrbitRegBaseExtended<MapT>>(
            VectorType{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
            { 0, -1, -1, 1,
              2,  3,  4,  5,
              -1, -1, -1, -1, 6 },
            { 0, 1, 2, 3, 7, 10 },

            std::cref(setup), t)
    {}
};

struct LyapunovOrbitRegParam
{
public:
    static RVector calculate(Real u0, size_t steps)
    {
        const Pcr3bp::SetupParameters<RMap> setup;

        const LyapunovOrbitRegLookupTable::Entry entry = calculate_initial(u0, steps, setup);

        LyapunovOrbitReg<RMap> orbit(setup, u0, g_t0);

        const RVector approx_root = { entry.pv0, entry.u1, entry.v1, entry.pu1, entry.pv1, entry.h1 };
        Carina::NewtonMethod newton(orbit, approx_root, steps);
        const RVector PV = newton.get_root();

        return RVector{ entry.u0, 0.0, 0.0, entry.pv0, entry.h1 };
    }

private:
    static LyapunovOrbitRegLookupTable::Entry calculate_initial(Real u0, size_t steps, const Pcr3bp::SetupParameters<RMap>& setup)
    {
        const LyapunovOrbitRegLookupTable::Entry base_entry
            = LyapunovOrbitRegLookupTable::get().get_base(u0);

        LyapunovOrbitReg_Implicit<RMap> orbit(setup, g_t0);

        Carina::ImplicitFunctionDifferentialEquation<decltype(orbit), 6> imp_fun(orbit);

        auto result = imp_fun.evolve( base_entry.u0,
            {
                Carina::scalar_cast<Real>(base_entry.pv0),
                Carina::scalar_cast<Real>(base_entry.u1),
                Carina::scalar_cast<Real>(base_entry.v1),
                Carina::scalar_cast<Real>(base_entry.pu1),
                Carina::scalar_cast<Real>(base_entry.pv1),
                Carina::scalar_cast<Real>(base_entry.h1)
            }, u0, steps);

        return LyapunovOrbitRegLookupTable::Entry{ u0, result[0], result[1], result[2], result[3], result[4], result[5] };
    }

    static constexpr Real g_t0 = 0.4;
};

}
