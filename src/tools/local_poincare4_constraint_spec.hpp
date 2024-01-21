///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/composite_map.hpp>
#include <capd_utils/local_map.hpp>
#include <capd_utils/gauss.hpp>

#include "local_poincare4_constraint_base.hpp"
#include "psi0_coefficients.hpp"

#include <proof/pcr3bp_reg_basic_objects.hpp>

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Implementation of the psi function (specialized for psi0)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LocalPoincare4_Constraint_Spec : public LocalPoincare4_Constraint_Base<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LocalPoincare4_Constraint_Spec(
        MapT& constraint,
        const CapdUtils::LocalCoordinateSystem<MapT>& src_coordsys)
            : LocalPoincare4_Constraint_Base<MapT>(constraint, src_coordsys)
            , m_src_coordsys_4_dim(src_coordsys)
    {}

    VectorType operator() (const VectorType& vec) override
    {
        return m_local_map(vec);
    }

    VectorType operator() (const VectorType& vec, MatrixType& der) override
    {
        return m_local_map(vec, der);
    }

    unsigned dimension() const override
    {
        return m_local_map.dimension();
    }

    unsigned imageDimension() const override
    {
        return m_local_map.imageDimension();
    }

private:
    const CapdUtils::LocalCoordinateSystem<MapT> m_src_coordsys_4_dim;

    MapT& m_internal_map
    {
        Psi0_Coefficients<MapT>::get().get_internal_map_ref()
    };

    static CapdUtils::LocalCoordinateSystem<MapT> create_src_coordsys(std::array<ScalarType, 2> d)
    {
        const ScalarType d1 = d.at(0);
        const ScalarType d2 = d.at(1);

        // print_var(d1);
        // print_var(d2);

        const VectorType origin = VectorType(2);
        MatrixType directions(2,2);
        directions(1,1) =  d1;
        directions(1,2) =  d1;
        directions(2,1) =  d2;
        directions(2,2) = -d2;

        return CapdUtils::LocalCoordinateSystem<MapT>(origin, directions);
    }

    const CapdUtils::LocalCoordinateSystem<MapT> m_src_coordsys_2_dim
    {
        create_src_coordsys( Psi0_Coefficients<MapT>::get().get_d_coeffs() )
    };
    
    CapdUtils::LocalMap<MapT,
        decltype(m_internal_map)&> m_local_map
    {
        std::ref(m_internal_map),
        std::ref(m_src_coordsys_2_dim),
        std::ref(m_src_coordsys_4_dim)
    };
};

}
