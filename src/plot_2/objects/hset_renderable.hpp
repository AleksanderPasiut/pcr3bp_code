///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "capd_renderable.hpp"

#include <capd_utils/local_coordinate_system.hpp>
#include <capd_utils/affine_map.hpp>
#include <capd_utils/enp_map.hpp>

#include <tools/local_poincare4_constraint.hpp>
#include <tools/local_poincare4_constraint_spec.hpp>
#include <tools/amortized_map.hpp>
#include <proof/pcr3bp_reg_basic_objects.hpp>

namespace Pcr3bpProof
{

class HsetRenderable
{
public:
    using MapT = RMap;
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    struct Params
    {
        Pcr3bp::RegBasicObjects<MapT>& basic_objects;
        CapdUtils::LocalCoordinateSystem<MapT> coordsys;
        std::array<double, 4> coordinates;
        size_t divs;
        size_t subdivs;
        float thickness;

        bool operator!= (const Params& params) const noexcept
        {
            return params.thickness != thickness;
        }
    };

    HsetRenderable(
        Lyra::Core3dObjects& core_objects_ref,
        const Manifold4_Transformation & transformation_ref,
        const Params& params)
            : m_params(params)
            , m_renderable(
                core_objects_ref,
                m_amortized_map,
                Leo::RulerSet<2>({
                    Leo::Ruler<double>( params.coordinates.at(0), params.coordinates.at(1), params.divs, params.subdivs ),
                    Leo::Ruler<double>( params.coordinates.at(2), params.coordinates.at(3), params.divs, params.subdivs )
                })
            , transformation_ref )
    {
        refresh();
    }

    void refresh()
    {
        m_renderable.fill(m_params.thickness);
    }

    const Params& get_params()
    {
        return m_params;
    }

    void show(bool arg)
    {
        m_renderable.show(arg);
    }

private:
    Params m_params;

    CapdUtils::AffineMap<MapT> m_linear
    {
        m_params.coordsys.get_origin(),
        m_params.coordsys.get_directions_matrix()
    };

    using LocalPoincare4_Constraint_BaseType = LocalPoincare4_Constraint_Base<MapT>;
    using LocalPoincare4_Constraint_BaseTypePtr = std::unique_ptr<LocalPoincare4_Constraint_BaseType>;
    using LocalPoincare4_Constraint_Type = LocalPoincare4_Constraint<MapT>;
    using LocalPoincare4_Constraint_SpecType = LocalPoincare4_Constraint_Spec<MapT>;

    LocalPoincare4_Constraint_BaseTypePtr m_extension_to_4_ptr
    {
        [this]() -> LocalPoincare4_Constraint_BaseTypePtr
        {
            auto origin = m_params.coordsys.get_origin();

            bool is_origin_zero = origin[0] == 0.0 && origin[1] == 0.0;

            return is_origin_zero ?
                LocalPoincare4_Constraint_BaseTypePtr(std::make_unique<LocalPoincare4_Constraint_SpecType>(
                    std::ref(m_params.basic_objects.m_hamiltonian_reg2),
                    std::ref(m_params.coordsys)
                )) :
                LocalPoincare4_Constraint_BaseTypePtr(std::make_unique<LocalPoincare4_Constraint_Type>(
                    std::ref(m_params.basic_objects.m_hamiltonian_reg2),
                    std::ref(m_params.coordsys)
                ));
        }()
    };

    LocalPoincare4_Constraint_BaseType& m_constraint
    {
        *m_extension_to_4_ptr
    };

    CapdUtils::CompositeMap<MapT, decltype(m_constraint)&, decltype(m_linear)&> m_composite
    {
        std::ref(m_constraint),
        std::ref(m_linear)
    };

    AmortizedMap<decltype(m_composite)> m_amortized_map
    {
        std::ref(m_composite)
    };

    CapdMapRenderable<decltype(m_amortized_map), 2, 4> m_renderable;
};

}
