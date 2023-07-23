///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"
#include "covering_relations_setup.hpp"

#include "tools/floating_info.hpp"

namespace Ursa
{

template<unsigned component_offset, typename LabelFunc>
void print_vector_list_tex(std::ostream& ostr, const std::list<IVector>& args, LabelFunc label_func)
{
    ostr << "\\renewcommand{\\arraystretch}{1.2}\n\n";
    ostr << "\\small\n\n";
    ostr << "\\begin{tabular}{ c | c | c }\n\n";
    ostr << "Vector & Component " << component_offset+1 << " & Component " << component_offset+2 << " \\\\\n\n";
    ostr << "\\hline\n\n";
    int index = 0;
    for (const IVector& vec : args)
    {
        ostr << "$";
        label_func(ostr, index);
        ostr << "$ & ";

        const Interval v0 = vec[component_offset];
        const Interval v1 = vec[component_offset+1];
        ostr << "$\\big[" << Carina::FloatingInfo<double>(v0.leftBound()) << ", " << Carina::FloatingInfo<double>(v0.rightBound()) << "\\big]$ & ";
        ostr << "$\\big[" << Carina::FloatingInfo<double>(v1.leftBound()) << ", " << Carina::FloatingInfo<double>(v1.rightBound()) << "\\big]$ \\\\\n\n";
        ++index;
    }
    ostr << "\\end{tabular}\n\n";
}

template<typename LabelFunc>
void print_vector_list_tex(std::ostream& ostr, const std::list<RVector>& args, LabelFunc label_func)
{
    ostr << "\\renewcommand{\\arraystretch}{1.2}\n\n";
    ostr << "\\small\n\n";
    ostr << "\\begin{tabular}{ c | c | c | c | c }\n\n";
    ostr << "Vector & Component 1 & Component 2 & Component 3 & Component 4 \\\\\n\n";
    ostr << "\\hline\n\n";
    int index = 0;
    for (const RVector& vec : args)
    {
        ostr << "$";
        label_func(ostr, index);
        ostr << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(vec[0]) << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(vec[1]) << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(vec[2]) << "$ & ";
        ostr << "$" << Carina::FloatingInfo<double>(vec[3]) << "$ \\\\\n\n";

        ++index;
    }
    ostr << "\\end{tabular}\n\n";
}

template<unsigned component_offset, typename LabelFunc>
void print_vector_list_tex(std::string filename, const std::list<IVector>& args, LabelFunc label_func)
{
    std::ofstream ostr(filename);

    if (ostr)
    {
        print_vector_list_tex<component_offset>(ostr, args, label_func);
        ostr.close();
    }
}

template<typename LabelFunc>
void print_vector_list_tex(std::string filename, const std::list<RVector>& args, LabelFunc label_func)
{
    std::ofstream ostr(filename);

    if (ostr)
    {
        print_vector_list_tex(ostr, args, label_func);
        ostr.close();
    }
}

TEST(Pcr3bp_proof, export_covering_relations_setup_data)
{
    using namespace Ursa;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};

    using Coordsys = Carina::LocalCoordinateSystem<IMap>;

    std::list<IVector> periodic_orbit_origins {};
    std::list<RVector> homoclinic_orbit_origins {};
    std::array<std::list<RVector>, 4> coordsys_directions {};

    for (Coordsys coordsys : setup.get_periodic_orbit_coordsys())
    {
        periodic_orbit_origins.emplace_back( coordsys.get_origin() );
        
        const std::array<IVector, 4> directions_iv
        {
            Carina::Extract<IMap>::get_vvector(coordsys.get_directions_matrix(), 1),
            Carina::Extract<IMap>::get_vvector(coordsys.get_directions_matrix(), 2),
            Carina::Extract<IMap>::get_vvector(coordsys.get_directions_matrix(), 3),
            Carina::Extract<IMap>::get_vvector(coordsys.get_directions_matrix(), 4)
        };

        for (size_t dir_idx = 0; dir_idx < directions_iv.size(); ++dir_idx)
        {
            const IVector direction_iv = directions_iv.at(dir_idx);
            assert_with_exception(Carina::span_vector(direction_iv) == RVector(4));

            coordsys_directions.at(dir_idx).emplace_back( Carina::vector_cast<RVector>(direction_iv) );
        }
    };

    for (Coordsys coordsys : setup.get_homoclinic_orbit_coordsys())
    {
        const IVector origin_iv = coordsys.get_origin();
        assert_with_exception(Carina::span_vector(origin_iv) == RVector(4));
        homoclinic_orbit_origins.emplace_back( Carina::vector_cast<RVector>(origin_iv) );

        const std::array<IVector, 4> directions_iv
        {
            Carina::Extract<IMap>::get_vvector(coordsys.get_directions_matrix(), 1),
            Carina::Extract<IMap>::get_vvector(coordsys.get_directions_matrix(), 2),
            Carina::Extract<IMap>::get_vvector(coordsys.get_directions_matrix(), 3),
            Carina::Extract<IMap>::get_vvector(coordsys.get_directions_matrix(), 4)
        };

        for (size_t dir_idx = 0; dir_idx < directions_iv.size(); ++dir_idx)
        {
            const IVector direction_iv = directions_iv.at(dir_idx);
            assert_with_exception(Carina::span_vector(direction_iv) == RVector(4));

            coordsys_directions.at(dir_idx).emplace_back( Carina::vector_cast<RVector>(direction_iv) );
        }
    };

    {
        auto label_func = [](std::ostream& ostr, unsigned index) { ostr << "c_{" << index << "}"; };
        print_vector_list_tex<0>("periodic_orbit_origins_1.tex.generated", periodic_orbit_origins, label_func);
    }
    {
        auto label_func = [](std::ostream& ostr, unsigned index) { ostr << "c_{" << index << "}"; };
        print_vector_list_tex<2>("periodic_orbit_origins_2.tex.generated", periodic_orbit_origins, label_func);
    }
    {
        auto label_func = [](std::ostream& ostr, unsigned index) { ostr << "c_{" << (index+4) << "}"; };
        print_vector_list_tex("homoclinic_orbit_origins.tex.generated", homoclinic_orbit_origins, label_func);
    }
    {
        auto label_func = [](std::ostream& ostr, unsigned index) { ostr << "w^{" << index << "}_1"; };
        print_vector_list_tex("coordsys_directions_1.tex.generated", coordsys_directions.at(0), label_func);
    }
    {
        auto label_func = [](std::ostream& ostr, unsigned index) { ostr << "w^{" << index << "}_2"; };
        print_vector_list_tex("coordsys_directions_2.tex.generated", coordsys_directions.at(1), label_func);
    }
    {
        auto label_func = [](std::ostream& ostr, unsigned index) { ostr << "w^{" << index << "}_3"; };
        print_vector_list_tex("coordsys_directions_3.tex.generated", coordsys_directions.at(2), label_func);
    }
    {
        auto label_func = [](std::ostream& ostr, unsigned index) { ostr << "w^{" << index << "}_4"; };
        print_vector_list_tex("coordsys_directions_4.tex.generated", coordsys_directions.at(3), label_func);
    }
}

}

