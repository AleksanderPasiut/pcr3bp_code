///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/capd/basic_types.hpp>
#include <capd_utils/map_base.hpp>

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Inverse coordinate change (from standard to regularized)
//! @details This coordinate change maps standard configuration space to fragment of regularized space for which u > 0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class LeviCivitaInverseCoordinateChange : public CapdUtils::MapBase<MapT>
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

	LeviCivitaInverseCoordinateChange(ScalarType xi, bool append_h = true)
        : m_xi(xi)
        , m_dimension(append_h ? 5 : 4)
        , m_f1(create_f1(xi, append_h))
        , m_f2(create_f2(xi, append_h))
        , m_f3(create_f3(xi, append_h))
        , m_f4(create_f4(xi, append_h))
        , m_g(create_g(append_h))
	{}

	VectorType operator()(const VectorType& x) override
	{
        this->assert_vector_size(x, m_dimension, "LeviCivitaInverseCoordinateChange x vector size mismatch (1)!");

		if(x[1] > 0)
        {
            return m_g( m_f1(x) );
        }
		if(x[1] < 0)
        {
            return m_g( m_f2(x) );
        }
		if(x[0] >= m_xi)
        {
            return m_f3(x);
        }
		return m_f4(x);
	}

    VectorType operator()(const VectorType& x, MatrixType& mat) override
	{
        this->assert_vector_size(x, m_dimension, "LeviCivitaInverseCoordinateChange x vector size mismatch (2)!");

		if(x[1] > 0)
        {
            MatrixType der1( m_f1.imageDimension(), m_f1.dimension() );
            const VectorType x1 = m_f1(x, der1);

            MatrixType der2 ( m_g.imageDimension(), m_g.dimension() );
            const VectorType x2 = m_g(x1, der2);

            mat = der2 * der1;
            return x2;
        }
		if(x[1] < 0)
        {
            MatrixType der1( m_f2.imageDimension(), m_f2.dimension() );
            const VectorType x1 = m_f2(x, der1);

            MatrixType der2 ( m_g.imageDimension(), m_g.dimension() );
            const VectorType x2 = m_g(x1, der2);

            mat = der2 * der1;
            return x2;
        }
		if(x[0] >= m_xi)
        {
            return m_f3(x, mat);
        }

		return m_f4(x, mat);
	}

    unsigned dimension() const noexcept override
    {
        return m_dimension;
    }

    unsigned imageDimension() const noexcept override
    {
        return m_dimension;
    }

private:
    static CapdUtils::Node delta(CapdUtils::Node& x, CapdUtils::Node& y, CapdUtils::Node& xi)
    {
        return sqrt( sqr(x-xi) + sqr(y) );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief y > 0 => u > 0 and v > 0
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MapT create_f1(ScalarType xi, bool append_h)
    {
        using CapdUtils::Node;
        auto func = [append_h](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& x = in[0];
            Node& y = in[1];
            Node& px = in[2];
            Node& py = in[3];
            Node& xi = param[0];

            out[0] = sqrt( ( delta(x, y, xi) + (x-xi) ) / 2 );
            out[1] = sqrt( ( delta(x, y, xi) - (x-xi) ) / 2 );
            out[2] = px;
            out[3] = py;
            
            if (append_h)
            {
                Node& h = in[4];
                out[4] = h;
            }
        };

        MapT map(func, m_dimension, m_dimension, 1);
        map.setParameter(0, xi);
        return map;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief y < 0 => u > 0 and v < 0
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MapT create_f2(ScalarType xi, bool append_h)
    {
        using CapdUtils::Node;
        auto func = [append_h](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& x = in[0];
            Node& y = in[1];
            Node& px = in[2];
            Node& py = in[3];
            Node& xi = param[0];

            out[0] =  sqrt( ( delta(x, y, xi) + (x-xi) ) / 2 );
            out[1] = -sqrt( ( delta(x, y, xi) - (x-xi) ) / 2 );
            out[2] = px;
            out[3] = py;
            
            if (append_h)
            {
                Node& h = in[4];
                out[4] = h;
            }
        };

        MapT map(func, m_dimension, m_dimension, 1);
        map.setParameter(0, xi);
        return map;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief y = 0 and x >= xi => u^2 = x-xi and v = 0
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MapT create_f3(ScalarType xi, bool append_h)
    {
        using CapdUtils::Node;
        auto func = [append_h](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& x = in[0];
            Node& y = in[1];
            Node& px = in[2];
            Node& py = in[3];

            Node& xi = param[0];

            out[0] = sqrt(x-xi);
            out[1] = Node(0);
            out[2] = 2*px*sqrt(x-xi);
            out[3] = 2*py*sqrt(x-xi);
            
            if (append_h)
            {
                Node& h = in[4];
                out[4] = h;
            }
        };

        MapT map(func, m_dimension, m_dimension, 1);
        map.setParameter(0, xi);
        return map;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief y = 0 and x < xi => u = 0 and v^2 = xi-x
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MapT create_f4(ScalarType xi, bool append_h)
    {
        using CapdUtils::Node;
        auto func = [append_h](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& x = in[0];
            Node& y = in[1];
            Node& px = in[2];
            Node& py = in[3];
            
            Node& xi = param[0];

            out[0] = Node(0);
            out[1] = sqrt(xi-x);
            out[2] =  2*py*sqrt(xi-x);
            out[3] = -2*px*sqrt(xi-x);
            
            if (append_h)
            {
                Node& h = in[4];
                out[4] = h;
            }
        };

        MapT map(func, m_dimension, m_dimension, 1);
        map.setParameter(0, xi);
        return map;
    };

    MapT create_g(bool append_h)
    {
        using CapdUtils::Node;
        auto func = [append_h](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& u = in[0];
            Node& v = in[1];
            Node& px = in[2];
            Node& py = in[3];

            out[0] = u;
            out[1] = v;
            out[2] = 2 * (px*u + py*v);
            out[3] = 2 * (py*u - px*v);

            if (append_h)
            {
                Node& h = in[4];
                out[4] = h;
            }
        };

        MapT map(func, m_dimension, m_dimension, 0);
        return map;
    }

    ScalarType m_xi;
    const unsigned m_dimension;

	MapT m_f1;
    MapT m_f2;
    MapT m_f3;
    MapT m_f4;
    MapT m_g;
};

}
