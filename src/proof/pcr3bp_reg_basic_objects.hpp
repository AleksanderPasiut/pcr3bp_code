///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/composite_map.hpp>

#include <pcr3bp_basic/setup_parameters.hpp>
#include <pcr3bp_basic/standard_system.hpp>
#include <pcr3bp_basic/regularized_system.hpp>

#include "periodic_orbit_parameters.hpp"

namespace Pcr3bpProof
{
namespace Pcr3bp
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief A container for the basic PCR3BP objects necessary for the computer-assisted proof
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class RegBasicObjects
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    Pcr3bp::SetupParameters<MapT> m_setup {};
    RegLyapunovCollisionOrbitParameters<MapT> m_parameters { m_setup };

    ScalarType m_h0 { m_parameters.get_energy() };

    MapT m_hamiltonian_reg2 { Pcr3bp::RegularizedSystem<MapT>::createHamiltonian4(2, m_setup, m_h0) };
    MapT m_hamiltonian_reg2_grad { Pcr3bp::RegularizedSystem<MapT>::createHamiltonianGradient4(2, m_setup, m_h0) };

    MapT m_vf_reg_pos2 { Pcr3bp::RegularizedSystem<MapT>::createPositiveVectorField4(2, m_setup, m_h0) };
    MapT m_vf_reg_neg2 { Pcr3bp::RegularizedSystem<MapT>::createNegativeVectorField4(2, m_setup, m_h0) };

    MapT m_collision_condition { Pcr3bp::RegularizedSystem<MapT>::createCollisionCondition(2, m_setup) };

    unsigned m_order { 60 };

    ScalarType m_lyapunov_orbit_period { 0.908942551524734 * 2 };

    struct ParallelogramCoveringsParameters
    {
        const ScalarType b0 { 255.0 / 256 }; // { 1.0 - pow(2.0, -24) };

        const ScalarType a0 { 151.0 / 256 };

        const ScalarType L { 0.000106 }; //{ 0.000105902 };

    } m_parallelogram_coverings_parameters;
};

}
}
