///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/map_base.hpp>

namespace Pcr3bpProof
{

template<typename MapT>
class LocalPoincare4_ProjectionBase : public CapdUtils::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    virtual VectorType operator() (const VectorType& vec) override = 0;

    virtual VectorType operator() (const VectorType& vec, MatrixType& der) override = 0;

    virtual unsigned dimension() const override = 0;

    virtual unsigned imageDimension() const override = 0;
};

}
