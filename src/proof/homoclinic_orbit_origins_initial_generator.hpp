///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author: Aleksander M. Pasiut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <capd_utils/enp_map.hpp>
#include <capd_utils/parallel_shooting/parallel_shooting.hpp>

#include <tools/auxiliary_functions.hpp>
#include <tools/local_poincare4_constraint.hpp>

#include "periodic_orbit_coordsys_generator.hpp"

namespace Pcr3bpProof
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief Generate initial homoclic orbit origins
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename MapT>
class HomoclinicOrbitOriginsInitialGenerator
{
public:
    using ScalarType = typename MapT::ScalarType;
    using VectorType = typename MapT::VectorType;
    using MatrixType = typename MapT::MatrixType;

    using Coordsys = CapdUtils::LocalCoordinateSystem<MapT>;

    HomoclinicOrbitOriginsInitialGenerator(const std::vector<Coordsys>& periodic_orbit_coordsys)
        : m_periodic_orbit_coordsys(periodic_orbit_coordsys)
    {
        const VectorType initial_root = m_epsmr_init( VectorType{ m_init_s } );
        CapdUtils::NewtonMethod newton( m_epsmr, initial_root, 100 );

        const VectorType root = newton.get_root();
        m_points = convert_root_into_initial_origins(root);

        m_total_expansion_factor = compute_total_expansion_factor_pos();
    }

    const std::vector<VectorType>& get_points() const noexcept
    {
        return m_points;
    }

    ScalarType get_total_expansion_factor() const noexcept
    {
        return m_total_expansion_factor;
    }

private:
    std::vector<VectorType> convert_root_into_initial_origins(VectorType root)
    {
        std::vector<VectorType> ret {};
        ret.reserve( 20 );

        for (int i = 1; i < root.dimension(); i += 3)
        {
            const VectorType v { root[i], 0.0, root[i+1], root[i+2] };
            ret.push_back(v);
        }

        const VectorType mid_v = m_poincare_pos_3( CapdUtils::Extract<MapT>::get_vector(root, root.dimension()-3, 3) );
        const VectorType mid_v4 = VectorType{ mid_v[0], 0.0, mid_v[1], mid_v[2] };
        ret.push_back(mid_v4);

        return ret;
    }

    ScalarType compute_total_expansion_factor_pos()
    {
        VectorType dir = CapdUtils::Extract<MapT>::get_vvector( m_coordsys_0.get_directions_matrix(), 1 );

        for (auto it = m_points.begin(); it != std::prev(m_points.end(), 1); ++it)
        {
            const VectorType origin_src = *it;

            MatrixType der(4,4);
            auto origin_img = m_poincare_pos(origin_src, der);

            dir = der * dir;
        }

        dir = AuxiliaryFunctions<MapT>::S_symmetry(dir);

        for (auto it = m_points.rbegin(); it != std::prev(m_points.rend(), 1); ++it)
        {
            const VectorType origin_src = *it;

            MatrixType der(4,4);
            auto origin_img = m_poincare_neg(origin_src, der);

            dir = der * dir;
        }

        return dir.euclNorm();
    }

    ScalarType compute_total_expansion_factor_neg()
    {
        VectorType dir = CapdUtils::Extract<MapT>::get_vvector( m_coordsys_0.get_directions_matrix(), 1 );

        for (auto it = m_points.rbegin(); it != m_points.rend(); ++it)
        {
            const VectorType origin_src = *it;

            MatrixType der(4,4);
            auto origin_img = m_poincare_neg(origin_src, der);

            dir = der * dir;
        }

        return dir.euclNorm();
    }

    const std::vector<Coordsys>& m_periodic_orbit_coordsys;

    Pcr3bp::RegBasicObjects<MapT> m_basic_objects {};

    CapdUtils::LocalCoordinateSystem<MapT> m_coordsys_0
    {
        m_periodic_orbit_coordsys.at(0)
    };

    CapdUtils::AffineMap<MapT> m_affine_0 { m_coordsys_0 };

    LocalPoincare4_Constraint<MapT> m_extension_to_4
    {
        m_basic_objects.m_hamiltonian_reg2,
        m_coordsys_0
    };

    CapdUtils::CompositeMap<MapT, MapT, decltype(m_extension_to_4)&, decltype(m_affine_0)&, MapT> m_init{ 
        CapdUtils::ExtensionMap<MapT>::create({ 0, -1 }),
        std::ref(m_extension_to_4),
        std::ref(m_affine_0),
        CapdUtils::ProjectionMap<MapT>::create(4, { 0, 2, 3 }) };


    CapdUtils::CoordinateSection<MapT> m_v_section { 4, 1, ScalarType(0.0) };
    CapdUtils::PoincareWrapper<MapT, decltype(m_v_section)> m_poincare_pos
    {
        m_basic_objects.m_vf_reg_pos2,
        m_basic_objects.m_order,
        m_v_section
    };

    CapdUtils::PoincareWrapper<MapT, decltype(m_v_section)> m_poincare_neg
    {
        m_basic_objects.m_vf_reg_neg2,
        m_basic_objects.m_order,
        m_v_section
    };

    CapdUtils::TimemapWrapper<MapT> m_timemap_pos
    {
        m_basic_objects.m_vf_reg_pos2,
        0.0,
        m_basic_objects.m_order
    };

    CapdUtils::ENP<MapT, decltype(m_poincare_pos)&> m_poincare_pos_3{ 
        { 0.0, 0.0, 0.0, 0.0 },
        { 0, -1, 1, 2 },
        { 0, 2, 3 },
        std::ref(m_poincare_pos)
    };

    CapdUtils::ENP<MapT, decltype(m_poincare_pos)&> m_poincare_pos_1{ 
        { 0.0, 0.0, 0.0, 0.0 },
        { 0, - 1, 1, 2 },
        { 2 },
        std::ref(m_poincare_pos)
    };

    ScalarType m_init_s {
        3.045e-9
    };

    CapdUtils::ParallelShootingInit<MapT, 
        decltype(m_init)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_1)&> m_epsmr_init {
            std::ref(m_init),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_1) };

    CapdUtils::ParallelShooting<MapT, 
        decltype(m_init)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_3)&,
        decltype(m_poincare_pos_1)&> m_epsmr {
            std::ref(m_init),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_3),
            std::ref(m_poincare_pos_1) };

    std::vector<VectorType> m_points {};

    ScalarType m_total_expansion_factor {};
};

}
