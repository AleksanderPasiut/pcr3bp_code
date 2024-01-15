///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/test_tools.hpp"
#include "covering_relations_setup.hpp"

TEST(Pcr3bp_proof, export_covering_relations_setup_data)
{
    using namespace Pcr3bpProof;

    capd::rounding::DoubleRounding::roundNearest();

    CoveringRelationsSetup setup {};

    using Coordsys = CapdUtils::LocalCoordinateSystem<IMap>;

    struct Data
    {
        RVector origin {};
        RVector unstable {};
        RVector stable {};
        RVector vector_field {};
        RVector energy {};

        Data(const Coordsys& coordsys)
        {
            const auto directions_matrix = CapdUtils::matrix_cast<RMatrix>(coordsys.get_directions_matrix());
            origin = CapdUtils::vector_cast<RVector>(coordsys.get_origin());
            unstable = CapdUtils::Extract<RMap>::get_vvector(directions_matrix, 1);
            stable = CapdUtils::Extract<RMap>::get_vvector(directions_matrix, 2);
            vector_field = CapdUtils::Extract<RMap>::get_vvector(directions_matrix, 3);
            energy = CapdUtils::Extract<RMap>::get_vvector(directions_matrix, 4);
        }
    };

    std::vector<Data> exported_data {};
    exported_data.reserve(18);

    for (const Coordsys& coordsys : setup.get_periodic_orbit_coordsys())
    {
        exported_data.emplace_back(std::cref(coordsys));
    }

    for (const Coordsys& coordsys : setup.get_homoclinic_orbit_coordsys())
    {
        exported_data.emplace_back(std::cref(coordsys));
    }

    auto get_origin = [](const Data& data) -> RVector { return data.origin; };
    auto get_unstable = [](const Data& data) -> RVector { return data.unstable; };
    auto get_stable = [](const Data& data) -> RVector { return data.stable; };
    auto get_vector_field = [](const Data& data) -> RVector { return data.vector_field; };
    auto get_energy = [](const Data& data) -> RVector { return data.energy; };

    auto export_latex_table = [&exported_data](
        const std::string& filename,
        std::function<RVector(const Data&)> getter,
        const std::string& label)
    {
        std::ofstream ofs(filename);
        ASSERT_TRUE(ofs);
        ofs.precision(6);

        int idx = 0;
        for (const Data& data : exported_data)
        {
            const RVector& v = getter(data);
            ofs << "$" << label << "_{" << idx << "}$ & $";
            ofs << v[0] << "$ & $";
            ofs << v[1] << "$ & $";
            ofs << v[2] << "$ & $";
            ofs << v[3] << "$ \\\\\n";
            ofs << "\\hline\n";
            ++idx;
        }

        ofs.close();
    };

    export_latex_table("output.origin.tex", get_origin, "w");
    export_latex_table("output.unstable.tex", get_unstable, "\\hat{u}");
    export_latex_table("output.stable.tex", get_stable, "\\hat{s}");
    export_latex_table("output.flow.tex", get_vector_field, "\\hat{v}");
    export_latex_table("output.energy.tex", get_energy, "\\hat{h}");
}
