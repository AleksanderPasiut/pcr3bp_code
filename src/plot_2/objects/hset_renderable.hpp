///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <lyra/core3d.hpp>
#include "capd_renderable.hpp"
#include "tools/test_tools.hpp"

#include <capd_utils/local_coordinate_system.hpp>
#include <capd_utils/affine_map.hpp>
#include <capd_utils/enp_map.hpp>

#include <tools/local_poincare4_constraint.hpp>
#include <tools/local_poincare4_constraint_spec.hpp>
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

    struct Param
    {
        Pcr3bp::RegBasicObjects<MapT>& basic_objects;
        CapdUtils::LocalCoordinateSystem<MapT> coordsys;
        std::array<double, 4> coordinates;
        size_t divs;
        size_t subdivs;
        const Manifold4_Transformation & transformation_ref;
        float thickness;
    };

    HsetRenderable(Lyra::Core3d& core_ref, const Param& param)
        : m_core_ref(core_ref)
        , m_param(param)
        , m_renderable(
            std::ref(core_ref.get_objects()),
            std::ref(m_composite),
            Leo::RulerSet<2>({
                Leo::Ruler<double>( param.coordinates.at(0), param.coordinates.at(1), param.divs, param.subdivs ),
                Leo::Ruler<double>( param.coordinates.at(2), param.coordinates.at(3), param.divs, param.subdivs ) }),
            std::cref(param.transformation_ref) )
    {
        refresh();
    }

    void refresh()
    {
        m_renderable.fill(m_param.thickness);
    }

private:
    Lyra::Core3d& m_core_ref;

    Param m_param;

    CapdUtils::AffineMap<MapT> m_linear
    {
        m_param.coordsys.get_origin(),
        m_param.coordsys.get_directions_matrix()
    };

    using LocalPoincare4_Constraint_BaseType = LocalPoincare4_Constraint_Base<MapT>;
    using LocalPoincare4_Constraint_BaseTypePtr = std::unique_ptr<LocalPoincare4_Constraint_BaseType>;
    using LocalPoincare4_Constraint_Type = LocalPoincare4_Constraint<MapT>;
    using LocalPoincare4_Constraint_SpecType = LocalPoincare4_Constraint_Spec<MapT>;

    LocalPoincare4_Constraint_BaseTypePtr m_extension_to_4_ptr
    {
        [this]() -> LocalPoincare4_Constraint_BaseTypePtr
        {
            auto origin = m_param.coordsys.get_origin();

            bool is_origin_zero = origin[0] == 0.0 && origin[1] == 0.0;

            return is_origin_zero ?
                LocalPoincare4_Constraint_BaseTypePtr(std::make_unique<LocalPoincare4_Constraint_SpecType>(
                    std::ref(m_param.basic_objects.m_hamiltonian_reg2),
                    std::ref(m_param.coordsys)
                )) :
                LocalPoincare4_Constraint_BaseTypePtr(std::make_unique<LocalPoincare4_Constraint_Type>(
                    std::ref(m_param.basic_objects.m_hamiltonian_reg2),
                    std::ref(m_param.coordsys)
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

    CapdMapRenderable<decltype(m_composite), 2, 4> m_renderable;
};

}
