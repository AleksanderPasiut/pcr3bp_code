///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <carina/capd/solution_curve.hpp>
#include <carina/capd/basic_tools.hpp>
#include <carina/concat.hpp>
#include <carina/gauss.hpp>
#include <carina/type_cast.hpp>

namespace Ursa
{

template<typename MapT>
class SolutionCurveWithCollisionCheck : public Carina::SolutionCurve<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using BoundType = typename ScalarType::BoundType;

    using CurvePieceType = typename Carina::SolutionCurve<MapT>::BaseCurve;

    SolutionCurveWithCollisionCheck() : Carina::SolutionCurve<MapT>(0.0)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @return True if condition is never satisfied. False if it might be satisfied.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool is_condition_never_satisfied(MapT condition, BoundType limit = 1e-15)
    {
        bool ret = true;
        
        for (CurvePieceType*& piece_ptr : this->pieces)
        {
            CurvePieceType& piece = *piece_ptr;

            ret &= internal_check(condition, limit, piece, piece.getLeftDomain(), piece.getRightDomain());
        }

        return ret;
    }

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Check if given condition is never satisfied for the specified curve piece
    //! @return True if condition is never satisfied. False if it might be satisfied.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static bool internal_check(MapT condition, BoundType limit, CurvePieceType& piece, BoundType left, BoundType right)
    {
        const ScalarType arg = ScalarType( left, right );
        const VectorType img = piece(arg);

        MatrixType condition_derivative( condition.imageDimension(), condition.dimension() );
        const VectorType val = condition(img, condition_derivative);

        if (image_does_not_intersect_with_zero(val))
        {
            return true;
        }

        if (Carina::span(arg) < limit)
        {
            print_var(arg);
            return false;
        }

        // Split the time interval and perform the check on parts recursively...
        const BoundType split_point = Carina::scalar_cast<BoundType>(arg);

        return
            internal_check( condition, limit, piece, left, split_point ) &&
            internal_check( condition, limit, piece, split_point, right );
    }

    static bool image_does_not_intersect_with_zero(VectorType image)
    {
        const VectorType zero_v = VectorType( image.dimension() );
        return capd::vectalg::intersectionIsEmpty( image, zero_v );
    }
};

}
