///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "setup_parameters.hpp"

namespace Pcr3bpProof
{
namespace Pcr3bp
{

template<typename MapT>
class StandardSystem
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Create Hamiltonian of PCR3BP in standard coordinates
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createHamiltonian(const Pcr3bp::SetupParameters<MapT>& setup, ScalarType h0 = ScalarType(0.0))
    {
        using CapdUtils::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& mu1 = param[0];
            Node& mu2 = param[1];
            Node& xi1 = param[2];
            Node& xi2 = param[3];
            Node& h0 = param[4];

            Node& x = in[0];
            Node& y = in[1];
            Node& px = in[2];
            Node& py = in[3];

            Node invr_1 = mu1*( (sqr(x-xi1) + sqr(y)) ^(-1.0/2) );
            Node invr_2 = mu2*( (sqr(x-xi2) + sqr(y)) ^(-1.0/2) );

            out[0] = ( sqr(px) + sqr(py) )/2 + y*px - x*py - invr_1 - invr_2 - h0;
        };
        
        MapT map(func, 4, 1, 5);
        map.setParameter(0, setup.get_mu(1));
        map.setParameter(1, setup.get_mu(2));
        map.setParameter(2, setup.get_x(1));
        map.setParameter(3, setup.get_x(2));
        map.setParameter(4, h0);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Create positive vector field of PCR3BP in standard coordinates
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createPositiveVectorField(const Pcr3bp::SetupParameters<MapT>& setup, bool h_coordinate = false)
	{
        return createVectorFieldInternal(setup, h_coordinate, ScalarType(+1.0));
	}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Create negative vector field of PCR3BP in standard coordinates
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createNegativeVectorField(const Pcr3bp::SetupParameters<MapT>& setup, bool h_coordinate = false)
	{
        return createVectorFieldInternal(setup, h_coordinate, ScalarType(-1.0));
	}

private:
    static MapT createVectorFieldInternal(const Pcr3bp::SetupParameters<MapT>& setup, bool h_coordinate, ScalarType direction)
	{
        using CapdUtils::Node;

        auto func = [h_coordinate](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& mu1 = param[0];
            Node& mu2 = param[1];
            Node& xi1 = param[2];
            Node& xi2 = param[3];
            Node& dir = param[4];

            Node& x = in[0];
            Node& y = in[1];
            Node& px = in[2];
            Node& py = in[3];

            out[0] = dir * (px + y);
            out[1] = dir * (py - x);

            Node invr_1 = mu1*( (sqr(x-xi1)+sqr(y))^(-3.0/2) );
		    Node invr_2 = mu2*( (sqr(x-xi2)+sqr(y))^(-3.0/2) );

            out[2] = dir * ( py - (x-xi1) * invr_1 - (x-xi2) * invr_2 );
            out[3] = dir * (-px - y * (invr_1 + invr_2) );

            if (h_coordinate)
            {
                Node& h = in[4];
                out[4] = Node(0.0);
            }
        };

        const int dimension = 4 + (h_coordinate ? 1 : 0);
        MapT map(func, dimension, dimension, 5);
        map.setParameter(0, setup.get_mu(1));
        map.setParameter(1, setup.get_mu(2));
        map.setParameter(2, setup.get_x(1));
        map.setParameter(3, setup.get_x(2));
        map.setParameter(4, direction);
		return map;
	}
};

}
}
