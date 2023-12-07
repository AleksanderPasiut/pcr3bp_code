///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/householder_generator2.hpp>
#include <capd_utils/local_coordinate_system.hpp>

#include "pcr3bp_reg_basic_objects.hpp"

namespace Pcr3bpProof
{
namespace Pcr3bp
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Generate local coordinate system in the specified point, where the basis directions (vu,vs,vf,vh) are orthonormal
//!        and:
//!        - vf is Pcr3bp reg2 vector field direction;
//!        - vh is Pcr3bp reg2 Hamiltonian gradient direction;
//!        - vu,vs are two unit vectors.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Reg2_InitialCoordsysGenerator
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    static CapdUtils::LocalCoordinateSystem<MapT> gen(Pcr3bp::RegBasicObjects<MapT>& basic_objects, VectorType point)
    {
        assert_with_exception(point.dimension() == 4);

        const VectorType vector_field_dir = basic_objects.m_vf_reg_pos2(point);
        const VectorType energy_dir = basic_objects.m_hamiltonian_reg2_grad(point);

        CapdUtils::HouseholderGenerator2<MapT> householder_generator2(vector_field_dir, energy_dir);

        const MatrixType dirs = householder_generator2.get_matrix();

        const MatrixType R
        {
            { 0, 0, 1, 0 },
            { 0, 0, 0, 1 },
            { 1, 0, 0, 0 },
            { 0, 1, 0, 0 }
        };

        return CapdUtils::LocalCoordinateSystem<MapT>(point, dirs * R);
    }
};

}
}
