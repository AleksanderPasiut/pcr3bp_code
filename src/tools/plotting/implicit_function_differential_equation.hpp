///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "capd/basic_types.hpp"
#include "capd/basic_tools.hpp"

#include "type_cast.hpp"

#include <leo/array/array_static_buffer.hpp>
#include <leo/diffeq.hpp>
#include <leo/linear/linear.hpp>

#include "map_base.hpp"

namespace Carina
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Implementation of implicit function differential equation
//!
//! @param MapT object implementing differential function itself;
//!             it must implement VectorType operator(VectorType, MatrixType&)
//!             which evaluates the map; the map must have size R^(N+1) -> R^N;
//! @param N    size of the implicit function argument vector
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT, size_t N>
class ImplicitFunctionDifferentialEquation
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using ArrayType = Leo::ArrayStaticBuffer<N, 1, Real>;

	ImplicitFunctionDifferentialEquation(MapT& map) : m_map(map)
	{
		if (m_map.imageDimension() != N)
		{
			throw std::invalid_argument("Map image dimension must be N!");
		}

		if (m_map.dimension() != N+1)
		{
			throw std::invalid_argument("Map dimension must be N+1!");
		}

        static_assert(Leo::is_read_indexable<ArrayType>());
	}

	VectorType evolve(ScalarType x0, VectorType p0, ScalarType x1, size_t steps, ScalarType dx_min = 1e-6) const
	{
		Leo::RK4_Step<ArrayType> rk4({N});

		ArrayType solution({N});

        for (size_t i = 0; i < N; ++i)
        {
            solution[{i}] = scalar_cast<Real>( p0[i] );
        }

		const Step step = get_step(x0, x1, steps, dx_min);
		Real x = scalar_cast<Real>(x0);
		for (size_t i = 0; i < step.count; ++i, x += step.dx)
		{
			rk4.move( solution, *this, x, step.dx);
		}

        VectorType ret(N);
        for (size_t i = 0; i < N; ++i)
        {
            ret[i] = solution[{i}];
        }

        return ret;
	}

private:
    struct Step
    {
        Real dx;
        size_t count;
    };

    static Step get_step(ScalarType x0, ScalarType x1, size_t count, ScalarType dx_min)
    {
        const Real delta_x = scalar_cast<Real>(x1 - x0);

        const size_t min_count = static_cast<size_t>(delta_x / scalar_cast<Real>(dx_min));

        Step ret {};
        ret.count = std::min(min_count+1, count);
        ret.dx = delta_x / ret.count;
        return ret;
    }

    template<typename ArrayLhs, typename ArrayRhs>
    void operator() (ArrayLhs& out, Real x, const ArrayRhs& arg) const
	{
        static_assert(Leo::is_write_indexable<ArrayLhs>());
        static_assert(Leo::is_read_indexable<ArrayLhs>());
        static_assert(Leo::is_read_indexable<ArrayRhs>());

        VectorType xp(N+1);
        xp[0] = x;

        for (size_t i = 0; i < N; ++i)
        {
            xp[i+1] = arg[{i}];
        }

        MatrixType der(N, N+1);
        m_map(xp, der);
        MapBase<MapT>::assert_matrix_size(der, N, N+1, "Implicit function matrix size mismatch!");

        if (N > 1)
        {
            Leo::ArrayStaticBuffer<N*N, 2, Real> mat({N, N});
            for (size_t i = 0; i < N; ++i)
            {
                for (size_t j = 0; j < N; ++j)
                {
                    mat[{i,j}] = scalar_cast<Real>( der(i+1, j+2) );
                }
            }

            for (size_t i = 0; i < N; ++i)
            {
                out[{i}] = -scalar_cast<Real>( der(i+1, 1) );
            }

            Leo::Linear::jordan_elimination(mat, out);
        }
        else
        {
    		out[{0}] = scalar_cast<Real>( - der(1, 1) / der(1, 2) );
        }
	}

    MapT& m_map;

    friend class Leo::RK4_Step<ArrayType>;
};

}
