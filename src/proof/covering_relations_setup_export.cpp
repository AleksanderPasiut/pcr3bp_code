///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"
#include "covering_relations_setup.hpp"

#include "tools/floating_info.hpp"

namespace Ursa
{

template<typename LabelFunc>
void print_vector_list_tex(std::ostream& ostr, const std::list<RVector>& args, LabelFunc label_func)
{
    ostr << "\\renewcommand{\\arraystretch}{1.1}\n\n";
    ostr << "\\small\n\n";
    ostr << "Vector & Component 1 & Component 2 & Component 3 & Component 4 \\\\\n\n";
    ostr << "\\begin{tabular}{ c | c | c | c | c }\n\n";
    ostr << "\\hline\n\n";
    int index = 0;
    for (const RVector& origin : args)
    {
        ostr << "$";
        label_func(ostr, index);
        ostr << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(origin[0]) << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(origin[1]) << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(origin[2]) << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(origin[3]) << "$ \\\\\n\n";

        ++index;
    }
    ostr << "\\end{tabular}\n\n";
}

}

TEST(Pcr3bp_proof, export_covering_relations_setup_data)
{
    using namespace Ursa;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};

    using Coordsys = Carina::LocalCoordinateSystem<IMap>;

    // std::ostream& ostr = std::cout;

    std::list<RVector> args {};
    for (Coordsys coordsys : setup.get_homoclinic_orbit_coordsys())
    {
        const IVector origin_iv = coordsys.get_origin();
        assert_with_exception(Carina::span_vector(origin_iv) == RVector(4));
        args.emplace_back( Carina::vector_cast<RVector>(origin_iv) );
    };

    auto label_func = [](std::ostream& ostr, unsigned index) { ostr << "c_{" << (index+4) << "}"; };

    std::ofstream ostr("homoclinic_orbit_origins.tex.generated");
    if (ostr)
    {
        print_vector_list_tex(ostr, args, label_func);
        ostr.close();
    }
    
    /*
    
    ostr << "\\renewcommand{\\arraystretch}{1.1}\n";
    ostr << "\\small\n";
    ostr << "Vector & Component 1 & Component 2 & Component 3 & Component 4 \\\\\n\n";
    ostr << "\\begin{tabular}{ c | c | c | c | c }\n\n";
    ostr << "\\hline\n\n";
    int index = 4;
    for (Coordsys coordsys : setup.get_homoclinic_orbit_coordsys())
    {
        const IVector origin_iv = coordsys.get_origin();
        assert_with_exception(Carina::span_vector(origin_iv) == RVector(4));

        const RVector origin = Carina::vector_cast<RVector>(origin_iv);

        ostr << "$c_{" << index << "}$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(origin[0]) << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(origin[1]) << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(origin[2]) << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(origin[3]) << "$ \\\\\n\n";

        ++index;
    }
    ostr << "\\end{tabular}\n\n";
    */

    // setup.get_periodic_orbit_coordsys()
}

