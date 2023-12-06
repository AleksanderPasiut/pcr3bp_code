///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "tools/test_tools.hpp"
#include "auxiliary_functions.hpp"

#include <capd_utils/local_coordinate_system.hpp>
#include <capd_utils/extract.hpp>
#include <capd_utils/concat.hpp>

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

    using Coordsys = CapdUtils::LocalCoordinateSystem<MapT>;

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

        const VectorType new_origin = { arg.get_origin()[0], 0.0, 0.0, arg.get_origin()[3] };

        const VectorType unstable_pos_4d = arg.get_directions_matrix() * VectorType{ unstable_pos_2d[0], unstable_pos_2d[1], 0.0, 0.0 };
        const VectorType unstable_neg_4d = AuxiliaryFunctions<MapT>::S_symmetry(unstable_pos_4d);

        const VectorType vf_dir = CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 3);
        const VectorType h_dir = CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 4);

        const VectorType vf_dir_s_symmetric = { 0.0, vf_dir[1], vf_dir[2], 0.0 };
        const VectorType h_dir_s_symmetric = { h_dir[0], 0.0, 0.0, h_dir[3] };

        const MatrixType new_dirs = CapdUtils::Concat<MapT>::build_matrix_from_vvectors({ unstable_pos_4d, unstable_neg_4d, vf_dir_s_symmetric, h_dir_s_symmetric });

        Coordsys ret
        {
            new_origin,
            new_dirs
        };

        return ret;
    }
};

}
