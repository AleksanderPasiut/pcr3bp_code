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
struct AdvancedConeConditions
{
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    AdvancedConeConditions(MatrixType mat)
    {
        u = Carina::scalar_cast<double>(mat(1,1));
        s = Carina::scalar_cast<double>(mat(2,2));

        auto abs = [](ScalarType iv) -> double
        {
            return std::max( std::abs(iv.leftBound()), std::abs(iv.rightBound()) );
        };

        c = std::max(
        {
            abs(mat(1,1) - u),
            abs(mat(1,2)),
            abs(mat(2,1)),
            abs(mat(2,2) - s)
        });

        EXPECT_TRUE(u>1.0);
        EXPECT_TRUE(s>0.0);
        EXPECT_TRUE(s<1.0);
        EXPECT_TRUE(c<s);
        EXPECT_TRUE(c>0.0);

        beta = (u - s)/c;
        
        p = (beta - sqrt(beta*beta - 4.0))/2 + 1e-6;

        EXPECT_TRUE(s > p*u);
        EXPECT_TRUE(p*(u-s) > c*(1+p)*(1+p));

        u_prim = (u - p*p*s) / (1 - p*p);
        s_prim = (s - p*u*u) / (1 - p*p);
        d = p*(u-s) / (1 - p*p);
        c_prim = c*(1+p)*(1+p) / (1 - p*p);


        EXPECT_TRUE( d + s_prim + 2*c_prim < 1.0);
        EXPECT_TRUE( s_prim + 4*c_prim < u_prim );

        alpha_min = (1 + d + c_prim) / (u_prim - c_prim);
        alpha_max = (1 - s_prim - c_prim) / (d + c_prim);

        alpha = alpha_min.rightBound() + 1e-6;

        std::ofstream fs("advanced_cone_condtion_parameters.txt");

        if (fs)
        {
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient u", u);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient s", s);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient c", c);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient beta", beta);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient p", p);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient u'", u_prim);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient s'", s_prim);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient c'", c_prim);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient d", d);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient alpha_min", alpha_min);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient alpha_max", alpha_max);
            Carina::VariablePrinter<MapT>::print(fs, "Coefficient alpha", alpha);

            fs.close();
        }
        else
        {
            throw std::logic_error("Failed to export parameters!");
        }
    }

    ScalarType u;
    ScalarType s;
    ScalarType c;
    ScalarType beta;

    ScalarType p;

    ScalarType u_prim;
    ScalarType s_prim;
    ScalarType d;
    ScalarType c_prim;

    ScalarType alpha_min;
    ScalarType alpha_max;
    ScalarType alpha;
};

}
