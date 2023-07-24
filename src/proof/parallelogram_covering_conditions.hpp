///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"
#include "tools/variable_printer.hpp"

#include <carina/type_cast.hpp>

namespace Ursa
{

template<typename MapT>
struct ParallelogramCoveringConditions
{
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    ParallelogramCoveringConditions(MatrixType mat)
    {
        phi = Carina::scalar_cast<double>(mat(1,1));
        theta = Carina::scalar_cast<double>(mat(2,2));

        auto abs = [](ScalarType iv) -> double
        {
            return std::max( std::abs(iv.leftBound()), std::abs(iv.rightBound()) );
        };

        c = std::max(
        {
            abs(mat(1,1) - phi),
            abs(mat(1,2)),
            abs(mat(2,1)),
            abs(mat(2,2) - theta)
        });

        EXPECT_TRUE(phi>1.0);
        EXPECT_TRUE(theta>0.0);
        EXPECT_TRUE(theta<1.0);
        EXPECT_TRUE(c<theta);
        EXPECT_TRUE(c>0.0);

        beta = (phi - theta)/c;
        
        p = (beta - sqrt(beta*beta - 4.0))/2 + 1e-6;

        EXPECT_TRUE(theta > p*phi);
        EXPECT_TRUE(p*(phi-theta) > c*(1+p)*(1+p));

        phi_prim = (phi - p*p*theta) / (1 - p*p);
        theta_prim = (theta - p*p*phi) / (1 - p*p);
        delta = p*(phi-theta) / (1 - p*p);
        c_prim = c*(1+p)*(1+p) / (1 - p*p);


        EXPECT_TRUE( delta + theta_prim + 2*c_prim < 1.0);
        EXPECT_TRUE( theta_prim + 4*c_prim < phi_prim );

        alpha_min = (theta_prim + delta + 2*c_prim) / (phi_prim - delta - 2*c_prim);
        alpha_max = (1 - theta_prim - c_prim) / (delta + c_prim);

        alpha = alpha_min.rightBound() + 1e-6;

        EXPECT_TRUE(alpha > 0.0);
        EXPECT_TRUE(alpha < 1.0);

        b_hat = alpha * (phi_prim - c_prim) - (delta+c_prim);

        std::ofstream fs("parallelogram_covering_conditions_parameters.txt");

        if (fs)
        {
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient phi", phi);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient theta", theta);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient c", c);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient beta", beta);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient p", p);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient phi'", phi_prim);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient theta'", theta_prim);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient c'", c_prim);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient delta", delta);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient alpha_min", alpha_min);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient alpha_max", alpha_max);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient alpha", alpha);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient b_hat", b_hat);

            fs.close();
        }
        else
        {
            throw std::logic_error("Failed to export parameters!");
        }
    }

    ScalarType phi;
    ScalarType theta;
    ScalarType c;
    ScalarType beta;

    ScalarType p;

    ScalarType phi_prim;
    ScalarType theta_prim;
    ScalarType delta;
    ScalarType c_prim;

    ScalarType alpha_min;
    ScalarType alpha_max;
    ScalarType alpha;

    ScalarType b_hat;
};

}
