///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "setup_parameters.hpp"

namespace Ursa
{
namespace Pcr3bp
{

template<typename MapT>
class RegularizedSystem
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Node = Carina::Node;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Create extended Hamiltonian of PCR3BP in regularized coordinates
    //!
    //! @param mu_index index of mass at which the regularization takes place
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createHamiltonian(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup)
    {
        auto func = [mu_index](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& mu_i = param[0];
            Node& mu3_i = param[1];
            Node& x_i = param[2];
            Node& epsilon = param[3];

            Node& u = in[0];
            Node& v = in[1];
            Node& pu = in[2];
            Node& pv = in[3];
            Node& h = in[4];

            Node factor = (sqr(sqr(u)-sqr(v)+epsilon) + sqr(2*u*v))^(-1.0/2);

            out[0] = Node(0);

            out[0] += 2*(sqr(u)+sqr(v))*(v*pu - u*pv - 2*mu3_i*factor - 2*h);
            out[0] -= 2*x_i*(v*pu+u*pv);
            out[0] -= 4*mu_i;
            out[0] += (sqr(pu) + sqr(pv)) / 2;
        };

        MapT map(func, 5, 1, 4);
        map.setParameter(0, setup.get_mu(mu_index));
        map.setParameter(1, setup.get_mu(3-mu_index));
        map.setParameter(2, setup.get_x(mu_index));
        map.setParameter(3, get_epsilon(mu_index));
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Create Hamiltonian of PCR3BP in regularized coordinates (fixed energy)
    //!
    //! @param mu_index index of mass at which the regularization takes place
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createHamiltonian4(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup, ScalarType h0)
    {
        auto func = [mu_index](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& mu_i = param[0];
            Node& mu3_i = param[1];
            Node& x_i = param[2];
            Node& epsilon = param[3];
            Node& h = param[4];

            Node& u = in[0];
            Node& v = in[1];
            Node& pu = in[2];
            Node& pv = in[3];

            Node factor = (sqr(sqr(u)-sqr(v)+epsilon) + sqr(2*u*v))^(-1.0/2);

            out[0] = Node(0);

            out[0] += 2*(sqr(u)+sqr(v))*(v*pu - u*pv - 2*mu3_i*factor - 2*h);
            out[0] -= 2*x_i*(v*pu+u*pv);
            out[0] -= 4*mu_i;
            out[0] += (sqr(pu) + sqr(pv)) / 2;
        };

        MapT map(func, 4, 1, 5);
        map.setParameter(0, setup.get_mu(mu_index));
        map.setParameter(1, setup.get_mu(3-mu_index));
        map.setParameter(2, setup.get_x(mu_index));
        map.setParameter(3, get_epsilon(mu_index));
        map.setParameter(4, h0);
        return map;
    }


    static ScalarType get_epsilon(size_t mu_index)
    {
        switch (mu_index)
        {
            case 1: return ScalarType(+1.0);
            case 2: return ScalarType(-1.0);
            default: throw std::logic_error("Unsupported i value!");
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Create positive vector field of PCR3BP in regularized coordinates
    //!
    //! @param mu_index index of mass at which the regularization takes place
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createPositiveVectorField(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup, bool t_coordinate)
    {
        return createVectorFieldInternal(mu_index, setup, ScalarType(+1.0), t_coordinate);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Create negative vector field of PCR3BP in regularized coordinates
    //!
    //! @param mu_index index of mass at which the regularization takes place
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createNegativeVectorField(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup, bool t_coordinate)
    {
        return createVectorFieldInternal(mu_index, setup, ScalarType(-1.0), t_coordinate);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Create gradient of Hamiltonian
    //!
    //! @param mu_index index of mass at which the regularization takes place
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createHamiltonianGradient4(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup, ScalarType h)
    {
        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& mu3_i = param[0];
            Node& x_i = param[1];
            Node& epsilon = param[2];
            Node& h = param[3];

            Node& u = in[0];
            Node& v = in[1];
            Node& pu = in[2];
            Node& pv = in[3];

            hamiltonianGradient(
                out[0],
                out[1],
                out[2],
                out[3],
                mu3_i,
                x_i,
                epsilon,
                h,
                u,
                v,
                pu,
                pv
                );
        };

        const int dimension = 4;
    
        MapT map(func, dimension, dimension, 4);
        map.setParameter(0, setup.get_mu(3-mu_index));
        map.setParameter(1, setup.get_x(mu_index));
        map.setParameter(2, get_epsilon(mu_index));
        map.setParameter(3, h);
        return map;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Create positive vector field of PCR3BP in regularized coordinates
    //!
    //! @param mu_index index of mass at which the regularization takes place
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createPositiveVectorField4(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup, ScalarType h)
    {
        return createVectorFieldInternal4(mu_index, setup, ScalarType(+1.0), h);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Create negative vector field of PCR3BP in regularized coordinates
    //!
    //! @param mu_index index of mass at which the regularization takes place
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static MapT createNegativeVectorField4(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup, ScalarType h)
    {
        return createVectorFieldInternal4(mu_index, setup, ScalarType(-1.0), h);
    }

private:
    static MapT createVectorFieldInternal(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup, ScalarType direction, bool t_coordinate)
    {
        auto func = [t_coordinate](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& mu3_i = param[0];
            Node& x_i = param[1];
            Node& epsilon = param[2];
            Node& dir = param[3];

            Node& u = in[0];
            Node& v = in[1];
            Node& pu = in[2];
            Node& pv = in[3];
            Node& h = in[4];

            hamiltonianGradient(
                out[2],
                out[3],
                out[0],
                out[1],
                mu3_i,
                x_i,
                epsilon,
                h,
                u,
                v,
                pu,
                pv
                );
            
            out[0] *= dir;
            out[1] *= dir;
            out[2] *= (-dir);
            out[3] *= (-dir);

            out[4] = Node(0);

            if (t_coordinate)
            {
                Node& t = in[5];
                out[5] = 4*(sqr(u)+sqr(v));
            }
        };

        const int dimension = 5 + (t_coordinate ? 1 : 0);
    
        MapT map(func, dimension, dimension, 4);
        map.setParameter(0, setup.get_mu(3-mu_index));
        map.setParameter(1, setup.get_x(mu_index));
        map.setParameter(2, get_epsilon(mu_index));
        map.setParameter(3, direction);
        return map;
    }

    static MapT createVectorFieldInternal4(size_t mu_index, const Pcr3bp::SetupParameters<MapT>& setup, ScalarType direction, ScalarType h)
    {
        auto func = [](Node, Node in[], int, Node out[], int, Node param[], int)
        {
            Node& mu3_i = param[0];
            Node& x_i = param[1];
            Node& epsilon = param[2];
            Node& dir = param[3];
            Node& h = param[4];

            Node& u = in[0];
            Node& v = in[1];
            Node& pu = in[2];
            Node& pv = in[3];

            hamiltonianGradient(
                out[2],
                out[3],
                out[0],
                out[1],
                mu3_i,
                x_i,
                epsilon,
                h,
                u,
                v,
                pu,
                pv
                );
            
            out[0] *= dir;
            out[1] *= dir;
            out[2] *= (-dir);
            out[3] *= (-dir);
        };

        MapT map(func, 4, 4, 5);
        map.setParameter(0, setup.get_mu(3-mu_index));
        map.setParameter(1, setup.get_x(mu_index));
        map.setParameter(2, get_epsilon(mu_index));
        map.setParameter(3, direction);
        map.setParameter(4, h);
        return map;
    }

    static void hamiltonianGradient(
        Node& ddu,
        Node& ddv,
        Node& ddpu,
        Node& ddpv,
        Node& mu3_i,
        Node& x_i,
        Node& epsilon,
        Node& h,
        Node& u,
        Node& v,
        Node& pu,
        Node& pv)
    {
        ddpu = pu + 2*v*(sqr(u) + sqr(v) - x_i);
        ddpv = pv - 2*u*(sqr(u) + sqr(v) + x_i);

        Node den = ( sqr(sqr(u)-sqr(v)+epsilon) + sqr(2*u*v) )^(-3.0/2);

        ddu = 4*u*(2*h-v*pu);
        ddu += 2*pv*(x_i+3*sqr(u)+sqr(v));
        ddu += 8*u*mu3_i*(1+epsilon*(sqr(u)-3*sqr(v))) * den;
        ddu *= -1;

        ddv = 4*v*(2*h+u*pv);
        ddv += 2*pu*(x_i-3*sqr(v)-sqr(u));
        ddv += 8*v*mu3_i*(1-epsilon*(sqr(v)-3*sqr(u))) * den;
        ddv *= -1;
    }
};

}
}
