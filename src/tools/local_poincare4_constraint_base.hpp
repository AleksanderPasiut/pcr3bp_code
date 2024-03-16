///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/map_base.hpp>
#include <capd_utils/local_coordinate_system.hpp>

#include "test_tools.hpp"

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Base class for psi function implementations
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LocalPoincare4_Constraint_Base : public CapdUtils::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    LocalPoincare4_Constraint_Base(
        MapT& constraint,
        const CapdUtils::LocalCoordinateSystem<MapT>& src_coordsys)
    {
        assert_with_exception(constraint.dimension() == 4);
        assert_with_exception(constraint.imageDimension() == 1);
        assert_with_exception(src_coordsys.get_origin().dimension() == 4);
        assert_with_exception(src_coordsys.get_origin().dimension() == 4);
    }

    virtual VectorType operator() (const VectorType& vec) override = 0;

    virtual VectorType operator() (const VectorType& vec, MatrixType& der) override = 0;

    virtual unsigned dimension() const override = 0;

    virtual unsigned imageDimension() const override = 0;
};

}
