///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/composite_map.hpp>
#include <capd_utils/local_map.hpp>
#include <capd_utils/gauss.hpp>

#include "local_poincare4_constraint_base.hpp"
#include "psi0_specialized.hpp"

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
            , m_dst_coordsys(src_coordsys)
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
    static CapdUtils::LocalCoordinateSystem<MapT> create_src_coordsys(
        MapT& internal_map,
        const CapdUtils::LocalCoordinateSystem<MapT>& dst_coordsys)
    {
        const VectorType unstable_dir = CapdUtils::Extract<MapT>::get_vvector(dst_coordsys.get_directions_matrix(), 1);
        print_var(unstable_dir);

        MatrixType dd(4, 2);
        internal_map( VectorType(2), dd );

        print_var(dd);

        MatrixType dd2(2,4);
        for (int i = 1; i <= 2; ++i)
            for (int j = 1; j <= 4; ++j)
                dd2(i,j) = dd(j,i);

        print_var(dd2);

        MatrixType ddx = dd2 * dd;

        print_var(ddx);

        const VectorType unstable_dir_reg = dd2 * unstable_dir;
        const VectorType d_coeff = CapdUtils::gauss<MapT>(ddx, unstable_dir_reg);

        const ScalarType d1 = d_coeff[0];
        const ScalarType d2 = d_coeff[1];

        const VectorType origin = VectorType(2);
        MatrixType directions(2,2);
        directions(1,1) =  d1;
        directions(1,2) =  d2;
        directions(2,1) =  d1;
        directions(2,2) = -d2;

        return CapdUtils::LocalCoordinateSystem<MapT>(origin, directions);
    }

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    MapT m_internal_map
    {
        Psi0_specialized<MapT>::create( m_basic_objects.m_h0, m_basic_objects.m_setup )
    };

    const CapdUtils::LocalCoordinateSystem<MapT> m_dst_coordsys;

    const CapdUtils::LocalCoordinateSystem<MapT> m_src_coordsys
    {
        create_src_coordsys(m_internal_map, m_dst_coordsys)
    };
    
    CapdUtils::LocalMap<MapT,
        decltype(m_internal_map)&> m_local_map
    {
        std::ref(m_internal_map),
        std::ref(m_src_coordsys),
        std::ref(m_dst_coordsys)
    };
};

}
