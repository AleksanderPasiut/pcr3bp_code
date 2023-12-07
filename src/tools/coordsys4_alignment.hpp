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
    //! @brief Replace first two vectors (columns) in the coordsys direction matrix with given vectors
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static Coordsys replace_unstable_dirs(Coordsys arg, VectorType unstable_pos, VectorType unstable_neg)
    {
        assert_with_exception(arg.get_origin().dimension() == 4);
        assert_with_exception(unstable_pos.dimension() == 4);
        assert_with_exception(unstable_neg.dimension() == 4);

        MatrixType new_directions_matrix = CapdUtils::Concat<MapT>::build_matrix_from_vvectors(
        {
            unstable_pos,
            unstable_neg,
            CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 3),
            CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 4)
        });

        Coordsys ret
        {
            arg.get_origin(),
            new_directions_matrix
        };

        return ret;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Replace first two vectors (columns) in the coordsys direction matrix with given vectors
    //! @details Ensures that the result coordinate system is s-symmetric
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static Coordsys replace_unstable_dirs_with_s_symmetry(Coordsys arg, VectorType unstable_pos)
    {
        assert_with_exception(arg.get_origin().dimension() == 4);
        assert_with_exception(unstable_pos.dimension() == 4);

        const VectorType new_origin = { arg.get_origin()[0], 0.0, 0.0, arg.get_origin()[3] };

        const VectorType unstable_neg = AuxiliaryFunctions<MapT>::S_symmetry(unstable_pos);

        const VectorType vf_dir = CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 3);
        const VectorType h_dir = CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 4);

        const VectorType vf_dir_s_symmetric = { 0.0, vf_dir[1], vf_dir[2], 0.0 };
        const VectorType h_dir_s_symmetric = { h_dir[0], 0.0, 0.0, h_dir[3] };

        MatrixType new_directions_matrix = CapdUtils::Concat<MapT>::build_matrix_from_vvectors(
        {
            unstable_pos,
            unstable_neg,
            vf_dir_s_symmetric,
            h_dir_s_symmetric
        });

        Coordsys ret
        {
            new_origin,
            new_directions_matrix
        };

        return ret;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! @brief Create coordinate system that is S-backsymmetric to the given one
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static Coordsys create_S_backsymmetric(Coordsys arg)
    {
        assert_with_exception(arg.get_origin().dimension() == 4);

        const VectorType u_dir = CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 1);
        const VectorType s_dir = CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 2);
        const VectorType v_dir = CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 3);
        const VectorType h_dir = CapdUtils::Extract<MapT>::get_vvector(arg.get_directions_matrix(), 4);

        MatrixType new_directions_matrix = CapdUtils::Concat<MapT>::build_matrix_from_vvectors(
        {
            AuxiliaryFunctions<MapT>::S_symmetry( s_dir ),
            AuxiliaryFunctions<MapT>::S_symmetry( u_dir ),
            -AuxiliaryFunctions<MapT>::S_symmetry( v_dir ),
            AuxiliaryFunctions<MapT>::S_symmetry( h_dir ),
        });

        Coordsys ret
        {
            AuxiliaryFunctions<MapT>::S_symmetry( arg.get_origin() ),
            new_directions_matrix
        };

        return ret;
    }
};

}
