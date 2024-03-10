///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <manifold4_transformation.hpp>
#include <pcr3bp_basic/setup_parameters.hpp>
#include <lyra/core3d.hpp>
#include <capd_utils/timemap_wrapper.hpp>

namespace Pcr3bpProof
{

class RegEvolution4 final
{
public:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    struct Params
    {
        Pcr3bp::SetupParameters<MapT> const & setup;
        VectorType initial_point;
        ScalarType t;
        size_t point_count;
        float thickness;
        bool positive;

        bool operator!= (const Params& arg) const noexcept
        {
            return
                this->initial_point != arg.initial_point ||
                this->t != arg.t ||
                this->point_count != arg.point_count ||
                this->thickness != arg.thickness;
        }
    };

    RegEvolution4(
        Lyra::Core3dObjects& core_objects_ref,
        const Manifold4_Transformation & transformation_ref,
        const Params& params);

    void refresh();

    const Params& get_params() const noexcept
    {
        return m_params;
    }

private:
    const Manifold4_Transformation & m_transformation_ref;

    MapT m_map;
    CapdUtils::TimemapWrapper<MapT> m_timemap;

    Params m_params;

    Leo::Ruler<ScalarType> m_ruler;

    Lyra::Manifold4<1> m_renderable;
};

}
