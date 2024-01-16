///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/capd/map.hpp>
#include <capd_utils/map_base.hpp>

#include <pcr3bp_basic/setup_parameters.hpp>


namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief A builder for the specialized psi0 function
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Psi0_specialized
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;
    
    static MapT create(const ScalarType& h0, const Pcr3bp::SetupParameters<MapT>& setup_parameters)
    {
        using CapdUtils::Node;

        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& x2 = param[0];
            Node& mu1 = param[1];
            Node& mu2 = param[2];
            Node& h = param[3];

            Node& u = in[0];
            Node& pu = in[1];

            out[0] = u;
            out[1] = Node(0.0);
            out[2] = pu;
            out[3] = 2*u*sqr(x2 + sqr(u)) + sqrt(
                4*sqr(u)*sqr(x2+sqr(u)) + 8*mu2 + 8*h*sqr(u) + 8*mu1*sqr(u) / sqrt(sqr(sqr(u) - 1)) - sqr(pu)
            );
        };

        MapT map(func, 2, 4, 4);
        map.setParameter(0, setup_parameters.get_x(2));
        map.setParameter(1, setup_parameters.get_mu(1));
        map.setParameter(2, setup_parameters.get_mu(2));
        map.setParameter(3, h0);
		return map;
    }
};

}
