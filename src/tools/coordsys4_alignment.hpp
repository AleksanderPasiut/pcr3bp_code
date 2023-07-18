///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"
#include "auxiliary_functions.hpp"

#include <carina/local_coordinate_system.hpp>
#include <carina/extract.hpp>
#include <carina/concat.hpp>

namespace Ursa
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief A utility that modifies the specified coordinate system according to the requirements
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class Coordsys4_Alignment
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = Carina::LocalCoordinateSystem<MapT>;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Align the system by replacing first two of its directions with unstable directions
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static Coordsys align(Coordsys arg, VectorType unstable_pos_2d, VectorType unstable_neg_2d)
    {
        assert_with_exception(arg.get_origin().dimension() == 4);
        assert_with_exception(unstable_pos_2d.dimension() == 2);
        assert_with_exception(unstable_neg_2d.dimension() == 2);

        MatrixType alignment_multiplier = MatrixType::Identity(4);
        alignment_multiplier(1,1) = unstable_pos_2d(1);
        alignment_multiplier(2,1) = unstable_pos_2d(2);
        alignment_multiplier(1,2) = unstable_neg_2d(1);
        alignment_multiplier(2,2) = unstable_neg_2d(2);

        Coordsys ret
        {
            arg.get_origin(),
            arg.get_directions_matrix() * alignment_multiplier
        };

        return ret;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Align the system by replacing first two of its directions with unstable directions
    //! @details Ensures that the result coordinate system is s-symmetric
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static Coordsys align_with_s_symmetry(Coordsys arg, VectorType unstable_pos_2d)
    {
        assert_with_exception(arg.get_origin().dimension() == 4);
        assert_with_exception(unstable_pos_2d.dimension() == 2);
        assert_with_exception(arg.get_origin()[1] == 0.0);
        assert_with_exception(arg.get_origin()[2] == 0.0);

        const VectorType unstable_pos_4d = arg.get_directions_matrix() * VectorType{ unstable_pos_2d[0], unstable_pos_2d[1], 0.0, 0.0 };
        const VectorType unstable_neg_4d = -AuxiliaryFunctions<MapT>::S_symmetry(unstable_pos_4d);

        // const VectorType vf_dir = Carina::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 3);
        const VectorType h_dir = Carina::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 4);

        // const VectorType vf_dir_s_symmetric = { vf_dir[0], 0.0, 0.0, vf_dir[3] };
        const VectorType h_dir_s_symmetric = { h_dir[0], 0.0, 0.0, h_dir[3] };

        const MatrixType new_dirs = Carina::Concat<MapT>::build_matrix_from_vvectors({ unstable_pos_4d, unstable_neg_4d, VectorType{ 0, 1, 0, 0 }, h_dir_s_symmetric });

        Coordsys ret
        {
            arg.get_origin(),
            new_dirs
        };

        return ret;
    }
};

}
