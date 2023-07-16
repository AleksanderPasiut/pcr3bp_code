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

    SolutionCurveWithCollisionCheck() : Carina::SolutionCurve<MapT>(0.0)
    {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @return True if condition is never satisfied. False if it might be satisfied.
    //! @details Implementation makes use of Newton theorem
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool is_condition_never_satisfied(MapT condition, VectorType x_arg, MatrixType e_der, BoundType limit = 1e-15)
    {
        bool ret = true;
        
        for (auto& piece_ptr : this->pieces)
        {
            auto& piece = *piece_ptr;

            auto internal_check = [&](double left, double right, auto& internal_check_ref) -> bool
            {
                const ScalarType arg = ScalarType( left, right );
                const VectorType img = piece(arg);

                MatrixType c_der( condition.imageDimension(), condition.dimension() );
                const VectorType val = condition(img, c_der);

                // print_var(arg);
                // print_var(val);
                const VectorType zero_v = VectorType( condition.imageDimension() );
                const bool intersection_empty = capd::vectalg::intersectionIsEmpty( val, zero_v );
                if (intersection_empty)
                {
                    return true;
                }

                try
                {
                    const VectorType t_dr = piece.timeDerivative(arg);
                    const MatrixType x_dr = piece.derivative(arg) * e_der;

                    MatrixType der = MatrixType( x_dr.dimension().first, x_dr.dimension().second+1 );
                    Carina::Concat<MapT>::copy_matrix_on_matrix(der, x_dr, 0, 1);

                    for (unsigned i = 1; i <= t_dr.dimension(); ++i)
                    {
                        der(i, 1) = t_dr(i);
                    }

                    // print_var(arg);
                    // print_var(x_arg);

                    const MatrixType m = c_der * der;
                    // print_var(m);

                    const MatrixType m_inv = Carina::gaussInverseMatrix<MapT>(m);
                    // print_var(m_inv);

                    const VectorType X0 = VectorType{ arg, x_arg[0], x_arg[1] };
                    const VectorType x0 = Carina::mid_vector(X0);

                    const VectorType newton = x0 - m_inv * Carina::mid_vector(val);
                    print_var(newton);

                    if (capd::vectalg::intersectionIsEmpty(X0, newton))
                    {
                        return true;
                    }

                    if (capd::vectalg::subset(newton, X0))
                    {
                        return true;
                    }
                }
                catch (const std::exception& e)
                {
                    // std::cout << e.what();
                }

                if (Carina::span(arg) < limit)
                {
                    print_var(arg);
                    return false;
                }

                const BoundType split = Carina::scalar_cast<BoundType>(arg);

                return 
                    internal_check_ref( left, split, internal_check_ref) &&
                    internal_check_ref( split, right, internal_check_ref);
            };

            ret &= internal_check(piece.getLeftDomain(), piece.getRightDomain(), internal_check);
        }

        return ret;
    }
};

}
