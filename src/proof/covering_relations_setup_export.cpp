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

    auto export_coordsys = [](std::ostream& ostr, const Coordsys& coordsys, int N)
    {
        ostr << "N" << N << ';';

        const auto origin = CapdUtils::vector_cast<RVector>(coordsys.get_origin());
        const auto directions_matrix = CapdUtils::matrix_cast<RMatrix>(coordsys.get_directions_matrix());
        for (int i = 0; i < 4; ++i)
        {
            ostr << origin[i] << ';';
        }

        for (int k = 0; k < 4; ++k)
        {
            const RVector direction = CapdUtils::Extract<RMap>::get_vvector(directions_matrix, k+1);

            for (int i = 0; i < 4; ++i)
            {
                ostr << direction[i] << ';';
            }
        }

        ostr << '\n';
    };

    std::ofstream ofs("output.csv");

    ASSERT_TRUE( bool(ofs) );

    ofs.precision(16);

    int index = 0;

    for (Coordsys coordsys : setup.get_periodic_orbit_coordsys())
    {
        export_coordsys(ofs, coordsys, index);
        ++index;
    }

    for (Coordsys coordsys : setup.get_homoclinic_orbit_coordsys())
    {
        export_coordsys(ofs, coordsys, index);
        ++index;
    }
}

